/ parser.q — standalone Barchart options-history CSV loader & normaliser
/ Namespace: .parser.barchart (public) / .parser.barchart.__ (private)
/ Standalone — no dependency on qFDM library modules
/ CamelCase column naming throughout

/ ══════════════════════════════════════════════════════════════════
/ CONFIGURATION
/ ══════════════════════════════════════════════════════════════════

.parser.barchart.cfg.fnPrefix:      "price-history-for-";
.parser.barchart.cfg.strikeDivisor: 100f;
.parser.barchart.cfg.metaCols:      `sym`expiry`strikeVal`cpFlag;
.parser.barchart.cfg.contractMult:  100f;

/ ══════════════════════════════════════════════════════════════════
/ PRIVATE HELPERS (__)
/ ══════════════════════════════════════════════════════════════════

/ Convert "YYYYMMDD" string to q date (for filename expiry parsing)
.parser.barchart.__toDate:{[s] "D"$(4#s),".",(s 4 5),".",s 6 7};

/ Convert Barchart Time column string to q date
/ Handles: "MM/DD/YY", "MM/DD/YYYY", "YYYY-MM-DD"
.parser.barchart.__parseTimeStr:{[s]
    if[10=count s;
       if[s[4]="-"; :"D"$s]];
    parts:"/" vs s;
    if[3=count parts;
       mm:parts 0; dd:parts 1; yy:parts 2;
       yyyy:$[2=count yy;"20",yy;yy];
       :"D"$yyyy,".",mm,".",dd];
    0Nd
 };

/ Clean one CSV header into a camelCase q symbol
/ "Open Int" -> `openInt  "%Change" -> `pctChange  "Price~" -> `price  "Time" -> `tradeDate
.parser.barchart.__cleanCol:{[rawName]
    cleaned:rawName except "\"~";
    cleaned:ssr[cleaned;"%";"Pct"];
    if[cleaned~"Time"; :`tradeDate];
    words:" " vs cleaned;
    firstWord:lower first words;
    restWords:{(upper 1#x),(lower 1_x)} each 1_words;
    `$firstWord,raze restWords
 };

/ Extract contract metadata from filename
/ price-history-for-aapl_20240216_16500c-05-27-2026.csv
/ -> sym:`aapl  expiry:2024.02.16  strikeVal:165.0  cpFlag:`c
.parser.barchart.__extractMeta:{[fpath]
    fpath:ssr[fpath;"\\";"/"];
    nm:first "." vs last "/" vs fpath;
    parts:"_" vs nm;
    sym:`$last "-" vs parts 0;
    expDt:.parser.barchart.__toDate parts 1;
    scp:first "-" vs parts 2;
    cpFlag:`$(-1)#scp;
    strikeRaw:(-1)_scp;
    strikeVal:("J"$strikeRaw)%.parser.barchart.cfg.strikeDivisor;
    .parser.barchart.cfg.metaCols!(sym;expDt;strikeVal;cpFlag)
 };

/ Column type sniffer — takes a string vector, returns typed vector
.parser.barchart.__sniffCol:{[v]
    nonEmpty:v where 0<count each v;
    if[0=count nonEmpty; :v];
    sample:first nonEmpty;
    / Date patterns
    if[sample like "[01][0-9]/[0-3][0-9]/[0-9]*"; :.parser.barchart.__parseTimeStr each v];
    if[sample like "[12][0-9][0-9][0-9]-[01][0-9]-[0-3]*"; :"D"$v];
    / Percent strings like "+1.23%" or "27.55%"
    if[any "%"=sample; :{"F"$x except "%+"} each v];
    / Try float parse (covers decimals, negatives, N/A -> null)
    fv:"F"$v;
    hasContent:0<count each v;
    validAsFloat:(not null fv) or not hasContent;
    if[all validAsFloat; :fv];
    / Fallback: symbols
    `$v
 };

/ Read one Barchart CSV file — returns typed table with metadata columns
.parser.barchart.__readOne:{[fpath]
    fpath:ssr[fpath;"\\";"/"];
    raw:read0 hsym `$fpath;
    / Remove footer rows, quoted lines, empty lines
    raw:raw where not raw like "\"*";
    raw:raw where not raw like "Downloaded*";
    raw:raw where 0<count each raw;
    if[2>count raw; '"empty or header-only file: ",fpath];
    / Parse header
    hdr:.parser.barchart.__cleanCol each "," vs first raw;
    / Parse body — keep only rows starting with date-like patterns
    bodyRows:1_raw;
    bodyRows:bodyRows where (bodyRows like "[01][0-9]/[0-3]*") or bodyRows like "[12][0-9][0-9][0-9]-*";
    if[0=count bodyRows; '"no valid data rows in: ",fpath];
    / Split and type columns
    splitRows:"," vs/: bodyRows;
    colVecs:flip splitRows;
    if[not (count hdr)=count colVecs; '"column count mismatch in: ",fpath];
    typedCols:.parser.barchart.__sniffCol each colVecs;
    tbl:flip hdr!typedCols;
    / Attach filename metadata
    md:.parser.barchart.__extractMeta fpath;
    nRows:count bodyRows;
    tbl,'flip .parser.barchart.cfg.metaCols!(nRows#'md .parser.barchart.cfg.metaCols)
 };

/ Safe wrapper for __readOne — returns empty list on error, logs warning
.parser.barchart.__safeReadOne:{[fpath]
    @[.parser.barchart.__readOne;fpath;{[fp;e] -1 "  SKIP ",fp," — ",e; ()}[fpath;]]
 };

/ ══════════════════════════════════════════════════════════════════
/ PUBLIC API
/ ══════════════════════════════════════════════════════════════════

/ List matching CSV files under a root directory (Windows or Linux)
.parser.barchart.listFiles:{[rootDir]
    rootDir:ssr[rootDir;"\\";"/"];
    winDir:ssr[rootDir;"/";"\\"];
    rawPaths:@[system;"dir /s /b \"",winDir,"\\*.csv\" 2>NUL";{()}];
    if[0=count rawPaths;
       rawPaths:@[system;"find \"",rootDir,"\" -name \"*.csv\" -type f 2>/dev/null";{()}]];
    if[0=count rawPaths; '"no CSVs found under: ",rootDir];
    rawPaths:{ssr[x except "\r";"\\";"/"] } each rawPaths;
    rawPaths where rawPaths like ("*",.parser.barchart.cfg.fnPrefix,"*")
 };

/ Load all matching Barchart CSV files from a directory tree
.parser.barchart.loadAll:{[rootDir]
    fps:.parser.barchart.listFiles rootDir;
    nFiles:count fps;
    -1 "bcParser: loading ",string[nFiles]," file(s)...";
    loaded:.parser.barchart.__safeReadOne each fps;
    nonEmpty:loaded where 0<count each loaded;
    allRows:raze nonEmpty;
    if[0=count allRows; '"no rows loaded from any file"];
    -1 "bcParser: ",string[count allRows]," rows loaded from ",string[count nonEmpty]," of ",string[nFiles]," files";
    `tradeDate xasc allRows
 };

/ Normalise raw Barchart table into qFDM/backtest-ready schema
.parser.barchart.normalise:{[optionTable]
    nRows:count optionTable;
    / Core fields from raw columns
    snapshotDates:optionTable`tradeDate;
    underlyings:optionTable`sym;
    expiryDates:optionTable`expiry;
    strikeVals:optionTable`strikeVal;
    cpFlags:optionTable`cpFlag;
    / Option type from cpFlag
    optionTypes:nRows#`put;
    optionTypes[where cpFlags=`c]:`call;
    / Contract ID (loop to avoid each-on-table gotcha)
    contractIds:nRows#`;
    cidIdx:0;
    while[cidIdx<nRows;
          contractIds[cidIdx]:`$string[underlyings cidIdx],"_",string[expiryDates cidIdx],"_",string[strikeVals cidIdx],"_",string cpFlags cidIdx;
          cidIdx+:1];
    / Spot (Barchart "Price~" column, cleaned to `price)
    spots:$[`price in cols optionTable;`float$optionTable`price;nRows#0Nf];
    / OHLC and quotes
    openVals:$[`open in cols optionTable;`float$optionTable`open;nRows#0Nf];
    highVals:$[`high in cols optionTable;`float$optionTable`high;nRows#0Nf];
    lowVals:$[`low in cols optionTable;`float$optionTable`low;nRows#0Nf];
    latestVals:$[`latest in cols optionTable;`float$optionTable`latest;nRows#0Nf];
    bidVals:$[`bid in cols optionTable;`float$optionTable`bid;nRows#0Nf];
    askVals:$[`ask in cols optionTable;`float$optionTable`ask;nRows#0Nf];
    / Mid and market price
    midVals:(bidVals+askVals)%2f;
    mktPrices:latestVals;
    nullLatest:where null latestVals;
    mktPrices[nullLatest]:midVals nullLatest;
    / IV: Barchart exports as percentage (27.55 = 27.55%), divide by 100
    rawIv:$[`iv in cols optionTable;`float$optionTable`iv;nRows#0Nf];
    impliedVols:rawIv%100f;
    / Volume and open interest
    volumeVals:$[`volume in cols optionTable;`float$optionTable`volume;nRows#0Nf];
    oiVals:$[`openInt in cols optionTable;`float$optionTable`openInt;nRows#0Nf];
    / Vendor Greeks (prefixed to distinguish from qFDM model Greeks)
    vDelta:$[`delta in cols optionTable;`float$optionTable`delta;nRows#0Nf];
    vGamma:$[`gamma in cols optionTable;`float$optionTable`gamma;nRows#0Nf];
    vTheta:$[`theta in cols optionTable;`float$optionTable`theta;nRows#0Nf];
    vVega:$[`vega in cols optionTable;`float$optionTable`vega;nRows#0Nf];
    vRho:$[`rho in cols optionTable;`float$optionTable`rho;nRows#0Nf];
    vTheo:$[`theo in cols optionTable;`float$optionTable`theo;nRows#0Nf];
    / Derived fields
    dteVals:expiryDates-snapshotDates;
    tauVals:dteVals%365f;
    moneynessVals:strikeVals%spots;
    spreadVals:askVals-bidVals;
    spreadPctVals:spreadVals%midVals;
    / Status and error message
    statusVals:nRows#`OK;
    errorMsgs:nRows#enlist "";
    / Validation: marketPrice null
    nullMktIdx:where null mktPrices;
    statusVals[nullMktIdx]:`error;
    errorMsgs[nullMktIdx]:(count nullMktIdx)#enlist "missing marketPrice";
    / Validation: spot non-positive (only flag rows still OK)
    badSpotIdx:where ((null spots) or spots<=0f) and statusVals=`OK;
    statusVals[badSpotIdx]:`error;
    errorMsgs[badSpotIdx]:(count badSpotIdx)#enlist "non-positive spot";
    / Validation: tau non-positive
    badTauIdx:where (dteVals<=0i) and statusVals=`OK;
    statusVals[badTauIdx]:`error;
    errorMsgs[badTauIdx]:(count badTauIdx)#enlist "non-positive tau";
    / Validation: IV non-positive
    badIvIdx:where ((null impliedVols) or impliedVols<=0f) and statusVals=`OK;
    statusVals[badIvIdx]:`error;
    errorMsgs[badIvIdx]:(count badIvIdx)#enlist "non-positive impliedVolatility";
    / Build normalised table
    flip `snapshotDate`underlying`contractId`expiryDate`optionType`strike`spot`open`high`low`latest`bid`ask`mid`marketPrice`impliedVolatility`volume`openInterest`vendorDelta`vendorGamma`vendorTheta`vendorVega`vendorRho`vendorTheo`dte`tau`moneyness`bidAskSpread`bidAskSpreadPct`status`errorMessage!(
        snapshotDates;underlyings;contractIds;expiryDates;optionTypes;strikeVals;spots;
        openVals;highVals;lowVals;latestVals;bidVals;askVals;midVals;mktPrices;impliedVols;
        volumeVals;oiVals;vDelta;vGamma;vTheta;vVega;vRho;vTheo;
        dteVals;tauVals;moneynessVals;spreadVals;spreadPctVals;statusVals;errorMsgs)
 };

/ Group normalised history by contract
.parser.barchart.byContract:{[optionHist]
    `underlying`expiryDate`strike`optionType xgroup optionHist
 };

/ Summary statistics for normalised option history
.parser.barchart.summary:{[optionHist]
    statusCol:optionHist`status;
    okCnt:sum statusCol=`OK;
    errCnt:(count optionHist)-okCnt;
    contractIds:distinct optionHist`contractId;
    `rowCount`contractCount`okRows`errorRows`minDate`maxDate`minExpiry`maxExpiry!(
        count optionHist;count contractIds;okCnt;errCnt;
        min optionHist`snapshotDate;max optionHist`snapshotDate;
        min optionHist`expiryDate;max optionHist`expiryDate)
 };

/ One-day market PnL replay
/ For each contract present on both dates, computes:
/   marketPnl = quantity * 100 * (exitPrice - entryPrice)
.parser.barchart.oneDayReplay:{[optionHist;entryDate;exitDate;qty]
    / Entry snapshot: OK rows on entryDate
    entryRows:optionHist where ((optionHist`snapshotDate)=entryDate) and (optionHist`status)=`OK;
    if[0=count entryRows; '"no OK rows on entryDate ",string entryDate];
    / Exit snapshot: rows on exitDate
    exitRows:optionHist where (optionHist`snapshotDate)=exitDate;
    if[0=count exitRows; '"no rows on exitDate ",string exitDate];
    / Build exit lookup: contractId -> index
    exitCids:exitRows`contractId;
    exitLookup:exitCids!til count exitRows;
    / Match entry contracts with exit contracts
    entryCids:entryRows`contractId;
    matchMask:entryCids in exitCids;
    matchedEntries:entryRows where matchMask;
    if[0=count matchedEntries; '"no matching contracts between entry and exit dates"];
    / Build replay rows using each on the matched entry list (list of dicts, safe for each)
    replayFn:{[exitRows;exitLookup;entryDate;exitDate;qty;eRow]
        cid:eRow`contractId;
        xIdx:exitLookup cid;
        xRow:exitRows xIdx;
        entryMktPx:eRow`marketPrice;
        exitMktPx:xRow`marketPrice;
        priceDiff:exitMktPx-entryMktPx;
        mktPnl:qty*.parser.barchart.cfg.contractMult*priceDiff;
        spotMv:xRow[`spot]-eRow`spot;
        ivMv:xRow[`impliedVolatility]-eRow`impliedVolatility;
        `contractId`underlying`optionType`strike`expiryDate`entryDate`exitDate`quantity`entryMarketPrice`exitMarketPrice`marketPnl`entrySpot`exitSpot`spotMove`entryIV`exitIV`ivMove`status`errorMessage!(
            cid;eRow`underlying;eRow`optionType;eRow`strike;eRow`expiryDate;
            entryDate;exitDate;qty;entryMktPx;exitMktPx;mktPnl;
            eRow`spot;xRow`spot;spotMv;eRow`impliedVolatility;xRow`impliedVolatility;ivMv;
            `OK;"")
    };
    boundFn:replayFn[exitRows;exitLookup;entryDate;exitDate;qty;];
    resultRows:();
    rIdx:0;
    while[rIdx<count matchedEntries;
          resultRows:resultRows,enlist boundFn matchedEntries rIdx;
          rIdx+:1];
    resultRows
 };

/ Replay summary statistics
.parser.barchart.replaySummary:{[replay]
    statusCol:replay`status;
    okCnt:sum statusCol=`OK;
    errCnt:(count replay)-okCnt;
    pnlCol:replay`marketPnl;
    `tradeCount`okRows`errorRows`totalMarketPnl`worstMarketPnl`bestMarketPnl!(
        count replay;okCnt;errCnt;sum pnlCol;min pnlCol;max pnlCol)
 };

/ ══════════════════════════════════════════════════════════════════
/ CRUDE FUTURES-CURVE PARSER  (.parser.crude.*)  (v0.50)
/ ══════════════════════════════════════════════════════════════════
/ Turns Barchart WTI daily settle CSVs into (tenor, price) forward-curve
/ snapshots for the v0.51 commodity calibration. Standalone like the equity
/ parser above (no qFDM library dependency).
/ -----------------------------------------------------------------
/ FILES   : data/barchart/CRUDE/Day_CRUDE_YYYYMM00.csv (one per delivery month;
/           YYYYMM is the delivery year-month, trailing "00" is Barchart's day
/           filler). Hour_* files are NOT used.
/ COLUMNS : Time,Open,High,Low,Latest,Volume  (header present, no footer).
/           Time is ISO8601 with a DST-shifting offset (e.g.
/           2018-12-27T06:00:00+0000); the date is the first 10 chars only.
/           Latest is the SETTLE/close (the curve price). Volume is a long.
/ EXPIRY  : derived from the data, NOT a hardcoded CME rule -> a contract's
/           expiry (last-trade date) is the MAX date in its own file, and its
/           firstDate is the MIN date. (Verified: the Jan-2020 file ends
/           2019-12-19, the actual last trade.)
/ TESTS   : driven by synthetic in-memory CSV TEXT only; the real CSVs are
/           gitignored. examples/load_crude_curve.q reads the real folder.

.parser.crude.cfg.fnPrefix: "Day_CRUDE_";
.parser.crude.cfg.monthCodes: "FGHJKMNQUVXZ";

/ ── private helpers (__) ──────────────────────────────────────────

/ Normalise the loadContract input into a list of CSV lines. Accepts either a
/ multi-line CSV TEXT string (contains a newline), a single-line PATH string
/ (no newline -> read from disk), or an already-split list of lines.
.parser.crude.__toLines:{[csvTextOrPath]
    if[10h=type csvTextOrPath;
       :$[any csvTextOrPath="\n";
            "\n" vs csvTextOrPath except "\r";
            .parser.crude.__readFileLines csvTextOrPath]];
    csvTextOrPath
 };

.parser.crude.__readFileLines:{[path]
    read0 hsym `$ssr[path;"\\";"/"]
 };

/ Locate a required column by (lower-cased) header name; controlled error if absent.
.parser.crude.__findCol:{[lowerHeader;colName]
    idx:lowerHeader?colName;
    if[idx=count lowerHeader; '"crude parser: missing required column '",colName,"'"];
    idx
 };

/ ── public API ────────────────────────────────────────────────────

/ CL futures month code for a 1-12 month (F G H J K M N Q U V X Z).
.parser.crude.monthCode:{[deliveryMonth]
    if[(deliveryMonth<1)|deliveryMonth>12; '"crude parser: month out of range"];
    .parser.crude.cfg.monthCodes deliveryMonth-1
 };

/ Parse "Day_CRUDE_YYYYMM00.csv" -> (deliveryYear; deliveryMonth) as longs.
.parser.crude.contractMonthFromFilename:{[filename]
    fpath:ssr[filename;"\\";"/"];
    baseName:first "." vs last "/" vs fpath;
    parts:"_" vs baseName;
    if[3>count parts; '"crude parser: expected Day_CRUDE_YYYYMM00.csv, got ",filename];
    ymToken:last parts;
    if[8>count ymToken; '"crude parser: bad YYYYMM00 token in ",filename];
    deliveryYear:"J"$4#ymToken;
    deliveryMonth:"J"$ymToken 4 5;
    if[(deliveryMonth<1)|deliveryMonth>12; '"crude parser: invalid delivery month in ",filename];
    (deliveryYear;deliveryMonth)
 };

/ Load one contract from CSV text (or a file path) -> typed OHLCV+settle table
/ sorted by date. Schema: date, open, high, low, settle, volume. The contract's
/ expiry/firstDate are derived downstream as max/min of the date column.
.parser.crude.loadContract:{[csvTextOrPath]
    lines:.parser.crude.__toLines csvTextOrPath;
    lines:lines where 0<count each lines;
    if[2>count lines; '"crude parser: empty or header-only CSV"];
    header:"," vs first lines;
    lowerHeader:lower each header;
    timeIdx:.parser.crude.__findCol[lowerHeader;"time"];
    openIdx:.parser.crude.__findCol[lowerHeader;"open"];
    highIdx:.parser.crude.__findCol[lowerHeader;"high"];
    lowIdx:.parser.crude.__findCol[lowerHeader;"low"];
    settleIdx:.parser.crude.__findCol[lowerHeader;"latest"];
    volumeIdx:.parser.crude.__findCol[lowerHeader;"volume"];
    bodyRows:"," vs/: 1_lines;
    nCols:count header;
    bodyRows:bodyRows where nCols=count each bodyRows;
    if[0=count bodyRows; '"crude parser: no valid data rows"];
    colVecs:flip bodyRows;
    dateVec:"D"$10#/:colVecs timeIdx;
    settleVec:"F"$colVecs settleIdx;
    tbl:flip `date`open`high`low`settle`volume!(
        dateVec;
        "F"$colVecs openIdx;
        "F"$colVecs highIdx;
        "F"$colVecs lowIdx;
        settleVec;
        "J"$colVecs volumeIdx);
    `date xasc tbl
 };

/ Expiry / first-trade dates of a loaded contract table (derived from the data).
.parser.crude.contractExpiry:{[contractTable] max contractTable`date};
.parser.crude.contractFirstDate:{[contractTable] min contractTable`date};

/ List Day_CRUDE_*.csv files under a directory (Hour_* excluded), Win or Linux.
.parser.crude.listFiles:{[rootDir]
    rootDir:ssr[rootDir;"\\";"/"];
    winDir:ssr[rootDir;"/";"\\"];
    rawPaths:@[system;"dir /s /b \"",winDir,"\\*.csv\" 2>NUL";{()}];
    if[0=count rawPaths;
       rawPaths:@[system;"find \"",rootDir,"\" -name \"*.csv\" -type f 2>/dev/null";{()}]];
    if[0=count rawPaths; '"crude parser: no CSVs found under ",rootDir];
    rawPaths:{ssr[x except "\r";"\\";"/"]} each rawPaths;
    baseNames:{last "/" vs x} each rawPaths;
    rawPaths where baseNames like (.parser.crude.cfg.fnPrefix,"*")
 };

/ Per-file loader tagged with contractYM/expiry/firstDate, isolated by try-catch
/ so one bad file never crashes the whole load.
.parser.crude.__loadOneTagged:{[fpath]
    @[.parser.crude.__loadOneTaggedRaw;fpath;{[fp;e] -1 "  SKIP ",fp," - ",e; ()}[fpath;]]
 };

.parser.crude.__loadOneTaggedRaw:{[fpath]
    deliveryYM:.parser.crude.contractMonthFromFilename fpath;
    contractYM:(100*deliveryYM 0)+deliveryYM 1;
    tbl:.parser.crude.loadContract fpath;
    expiry:.parser.crude.contractExpiry tbl;
    firstDate:.parser.crude.contractFirstDate tbl;
    tagged:update contractYM:contractYM, expiry:expiry, firstDate:firstDate from tbl;
    `contractYM`expiry`firstDate`date`settle`open`high`low`volume xcols tagged
 };

/ Load all Day_CRUDE_*.csv under dir into one long table across contracts.
/ Schema: contractYM, expiry, firstDate, date, settle, open, high, low, volume.
.parser.crude.loadAll:{[rootDir]
    fps:.parser.crude.listFiles rootDir;
    -1 "crudeParser: loading ",string[count fps]," contract file(s)...";
    pieces:.parser.crude.__loadOneTagged each fps;
    nonEmpty:pieces where 0<count each pieces;
    if[0=count nonEmpty; '"crude parser: no rows loaded from any file"];
    longTable:raze nonEmpty;
    -1 "crudeParser: ",string[count longTable]," rows from ",string[count nonEmpty]," of ",string[count fps]," file(s)";
    `contractYM`date xasc longTable
 };

/ Forward curve as-of a date: (tenor, price) for every contract alive AND quoted
/ on asof (firstDate<=asof<=expiry and a settle row exists on asof). tenor =
/ (expiry-asof)/365 in years. Sorted ascending by tenor. This is the marketCurve
/ shape the v0.51 calibration consumes. Controlled error if nothing is alive.
.parser.crude.curveAt:{[longTable;asofArg]
    if[not 98h=type longTable; '"crude curveAt: expected a table"];
    asofDate:$[-14h=type asofArg; asofArg; "D"$asofArg];
    alive:select from longTable where firstDate<=asofDate, expiry>=asofDate, date=asofDate;
    if[0=count alive; '"crude curveAt: no contracts alive/quoted on ",string asofDate];
    grouped:0!select expiry:first expiry, price:first settle by contractYM from alive;
    grouped:update tenor:(`float$expiry-asofDate)%365f from grouped;
    `tenor xasc `tenor`price`contractYM`expiry xcols grouped
 };

.parser.crude.__curveAtSafe:{[longTable;asofArg]
    builder:{[lt;a] curve:.parser.crude.curveAt[lt;a]; update asofDate:a from curve};
    @[builder[longTable;];asofArg;{[e] ()}]
 };

/ Curve per asof date (curve evolution over time for calibration). Dates with no
/ alive contracts are skipped (not fatal). Schema: asofDate, tenor, price,
/ contractYM, expiry. (Column is asofDate, not asof: asof is a q built-in.)
.parser.crude.curveHistory:{[longTable;asofDates]
    asofVec:$[-14h=type asofDates; enlist asofDates; asofDates];
    pieces:.parser.crude.__curveAtSafe[longTable;] each asofVec;
    nonEmpty:pieces where 0<count each pieces;
    if[0=count nonEmpty; '"crude curveHistory: no curves on any asof date"];
    `asofDate`tenor`price`contractYM`expiry xcols raze nonEmpty
 };

-1 "parser.q loaded - .parser.barchart + .parser.crude namespaces ready";

/ ══════════════════════════════════════════════════════════════════
/ USAGE
/ ══════════════════════════════════════════════════════════════════
/ \l parser.q
/ optionRawTable:.parser.barchart.loadAll "C:/q/w64/qPricer/data/barchart/aapl/options_history";
/ optionTable:.parser.barchart.normalise optionTable;
/ .parser.barchart.summary optionHist
/ replay:.parser.barchart.oneDayReplay[optionHist;2024.01.02;2024.01.03;1f];
/ .parser.barchart.replaySummary replay

