/ data/hdb.q - partitioned HDB ingestion + query layer (.data.hdb.*) (v0.59)
/ ----------------------------------------------------------------------------
/ Migration step 3 (ARCHITECTURE.md sections 3 / 9). Replaces per-run CSV
/ re-parsing with a date-partitioned, splayed kdb+ HDB built once from the
/ Barchart CSVs. ADDITIVE: the parser (.parser.futures / .parser.crude) is
/ UNCHANGED and is now the CSV-read stage that FEEDS ingestion. The HDB only
/ replaces loadAll (the parse); the curve DERIVATION is REUSED verbatim from
/ .parser.crude.curveAt / curveHistory, so query output is byte-identical to the
/ parser's by construction (same long-table schema in -> same parser code).
/ -
/ Layout: ONE HDB for all commodities; table `futures` SPLAYED and UNPARTITIONED
/ (the dataset is ~48k rows total - tiny for kdb+ - so ~1900 daily partition dirs
/ would be slow and Windows-hostile; a single splayed dir is far faster). Columns:
/ commodity(sym, p# parted), contractYM(long, as the parser), expiry(date),
/ firstDate(date), date(date, a normal column - not a partition domain), open/
/ high/low/settle(float), volume(long). expiry/firstDate are the parser's
/ per-contract max/min dates denormalised onto each row, so the curve filter +
/ tenor arithmetic match exactly. Symbol columns are enumerated via .Q.en against
/ the sym file before the column files are written with `set`; queries are plain
/ `select from futures where commodity=c, date=d` (no partition domain).
/ -
/ Library-load independence: this file only DEFINES functions; it never opens an
/ HDB at import. The real HDB is gitignored and may be absent (fresh checkout /
/ CI / the synthetic test) - opening is explicit via .data.hdb.open.
/ ----------------------------------------------------------------------------

/ ── ingestion (touches real gitignored data; driven by scripts/ingest_hdb.q) ──

/ Remove the HDB directory so a rebuild is clean / idempotent. Branches on the OS
/ (.z.o) so only the platform-appropriate command runs - calling the POSIX form on
/ Windows cmd would mis-parse the flag as a path. Errors are swallowed (absent dir).
.data.hdb.__clean:{[hdbPath]
    $[.z.o like "w*";
        @[system;"rmdir /s /q \"",(ssr[hdbPath;"/";"\\"]),"\"";{[e]}];
        @[system;"rm -rf \"",hdbPath,"\"";{[e]}]];
    hdbPath
 };

/ Build the SPLAYED, UNPARTITIONED HDB at hdbPath from CSVs under csvDir for each
/ commodity tag in `commodities`. REUSES .parser.futures.loadAll verbatim (so
/ expiry = MAX date / firstDate = MIN date are derived exactly as the parser
/ does), tags each long table with its commodity, sorts by commodity (then
/ contractYM, date) and applies the p# attribute on commodity, enumerates sym
/ columns via .Q.en against hdbPath/sym, and writes the column files to a single
/ hdbPath/futures/ splay with `set`. Idempotent (the target is wiped first). No
/ mkdir is needed - .Q.en creates hdbPath + writes hdbPath/sym, and `set` creates
/ hdbPath/futures/; q creates intermediate directories on write.
.data.hdb.ingest:{[csvDir;commodities;hdbPath]
    commodities:(),commodities;
    .data.hdb.__clean hdbPath;
    hp:hsym `$hdbPath;
    longs:raze {[csvDir;tag]
        lt:.parser.futures.loadAll[csvDir;tag];
        update commodity:tag from lt
        }[csvDir;] each commodities;
    t:`commodity`contractYM`date xasc select commodity,contractYM,expiry,firstDate,date,open,high,low,settle,volume from longs;
    t:.Q.en[hp;t];
    t:update `p#commodity from t;
    (hsym `$hdbPath,"/futures/") set t;
    .Q.gc[];
    -1 "hdb.ingest: ",(string count t)," rows, ",(string count commodities),
       " commodit(ies), splayed -> ",hdbPath,"/futures/";
    `rows`commodities`dates!(count t;commodities;count distinct t`date)
 };

/ ── query layer (mirrors the parser's curve output exactly) ──

/ Open the splayed HDB: load the sym enumeration domain + the `futures` table
/ into the workspace. Uses `get` (NOT \l) deliberately - \l of a database
/ directory changes the process working directory to the HDB root, which would
/ break any later relative-path file access in a long-running process (e.g. the
/ test runner reading the next test file). `get` reads the splayed columns into
/ memory (trivial at this dataset size) with no cwd change and no lingering file
/ lock. Errors cleanly if the path / sym file is absent (the HDB is gitignored
/ and may not exist in a fresh checkout / CI).
.data.hdb.open:{[hdbPath]
    hp:hsym `$hdbPath;
    if[0=count key hp; '"hdb.open: path not found or empty: ",hdbPath];
    symPath:hsym `$hdbPath,"/sym";
    if[0=count key symPath; '"hdb.open: not a valid HDB (no sym file): ",hdbPath];
    sym::get symPath;
    futures::get hsym `$hdbPath,"/futures/";
    hdbPath
 };

.data.hdb.__requireLoaded:{[]
    if[not `futures in key `.; '"hdb: not opened - call .data.hdb.open[hdbPath] first"];
 };

/ Sorted distinct trade dates available for a commodity (mirrors
/ `asc distinct longTable`date` on the parser long table).
.data.hdb.dates:{[commodityArg]
    .data.hdb.__requireLoaded[];
    asc distinct exec date from futures where commodity=commodityArg
 };

/ Materialise the parser-schema long table for one commodity from the HDB
/ (one columnar read). Column names/types/order match .parser.futures.loadAll
/ output exactly, so it can be fed straight into the parser's curve functions.
.data.hdb.__longTable:{[commodityArg]
    .data.hdb.__requireLoaded[];
    select contractYM,expiry,firstDate,date,settle,open,high,low,volume from futures where commodity=commodityArg
 };

/ Forward curve as-of a date for one commodity. Filters the splayed table on
/ date + commodity (commodity uses the p# attribute), then reuses
/ .parser.crude.curveAt on the materialised rows -> (tenor, price, contractYM,
/ expiry) sorted by tenor, byte-identical to the parser.
.data.hdb.curveAt:{[commodityArg;asofArg]
    .data.hdb.__requireLoaded[];
    asofDate:$[-14h=type asofArg; asofArg; "D"$asofArg];
    lt:select contractYM,expiry,firstDate,date,settle,open,high,low,volume from futures where commodity=commodityArg, date=asofDate;
    .parser.crude.curveAt[lt;asofDate]
 };

/ Curve per asof date for one commodity (the panel calibration / Kalman /
/ backtest consume). One columnar read of the commodity, then the parser's
/ curveHistory verbatim -> byte-identical (asofDate, tenor, price, contractYM, expiry).
.data.hdb.curveHistory:{[commodityArg;asofDates]
    .parser.crude.curveHistory[.data.hdb.__longTable commodityArg;asofDates]
 };

-1 "hdb.q loaded - .data.hdb.* splayed HDB ingestion + query ready (not opened)";
