\l core/init.q
/ HDB ingest<->query equivalence on SYNTHETIC data (no real CSV/HDB). Writes
/ Barchart-format CSV text for two synthetic commodities under scratch/, ingests
/ into a temp SPLAYED HDB, and asserts .data.hdb.curveAt / curveHistory are
/ byte-identical (schema + values, incl. tenor arithmetic + expiry derivation) to
/ the parser's curveAt / curveHistory on the same input. Also checks the commodity
/ filter isolates correctly and that .data.hdb.open errors cleanly on a missing
/ path. This is the core correctness proof for migration step 3 and needs no real
/ data. NO cleanup: all scratch lives under scratch/ (gitignored) and is left in
/ place; leaving the mapped `futures` / `sym` globals is safe (no other test
/ references a bare `futures` global, and `sym` is only ever a local elsewhere).

mkdirRobust:{[d] $[.z.o like "w*"; @[system;"mkdir \"",(ssr[d;"/";"\\"]),"\"";{[e]}]; @[system;"mkdir -p \"",d,"\"";{[e]}]];};
hdr:"Time,Open,High,Low,Latest,Volume";
mkcsv:{[dts;settles] hdr,"\n","\n" sv {[d;s] (string d),"T06:00:00+0000,",(string s),",",(string s),",",(string s),",",(string s),",1000"}'[dts;settles]};

base:"scratch/hdb_eqv"; csvDir:base,"/csv"; hdbDir:base,"/hdb";
mkdirRobust csvDir;

dShort:2019.12.30 2019.12.31 2020.01.02 2020.01.03;
dLong:2019.12.30 2019.12.31 2020.01.02 2020.01.03 2020.01.06;
/ commodity TST: 3 contracts (one carries an extra later date to prove curveAt
/ picks the asof settle); commodity TST2: 2 contracts, different price level.
(hsym `$csvDir,"/Day_TST_20200100.csv") 0: "\n" vs mkcsv[dShort;50.0 50.1 50.2 50.3];
(hsym `$csvDir,"/Day_TST_20200200.csv") 0: "\n" vs mkcsv[dLong;51.0 51.1 51.2 51.3 51.4];
(hsym `$csvDir,"/Day_TST_20200300.csv") 0: "\n" vs mkcsv[dLong;52.0 52.1 52.2 52.3 52.4];
(hsym `$csvDir,"/Day_TST2_20200100.csv") 0: "\n" vs mkcsv[dShort;70.0 70.1 70.2 70.3];
(hsym `$csvDir,"/Day_TST2_20200200.csv") 0: "\n" vs mkcsv[dLong;71.0 71.1 71.2 71.3 71.4];

asofD:2020.01.02;
/ parser reference: per-commodity loadAll then parser curve funcs.
pltTst:.parser.futures.loadAll[csvDir;`TST];
pltTst2:.parser.futures.loadAll[csvDir;`TST2];
pCurveAtTst:.parser.futures.curveAt[pltTst;asofD];
pHistTst:.parser.futures.curveHistory[pltTst;asc distinct pltTst`date];
pCurveAtTst2:.parser.futures.curveAt[pltTst2;asofD];

/ HDB: ingest both commodities, open, query each.
ingRes:.data.hdb.ingest[csvDir;`TST`TST2;hdbDir];
.data.hdb.open hdbDir;
hCurveAtTst:.data.hdb.curveAt[`TST;asofD];
hHistTst:.data.hdb.curveHistory[`TST;.data.hdb.dates `TST];
hCurveAtTst2:.data.hdb.curveAt[`TST2;asofD];
hDatesTst:.data.hdb.dates `TST;

/ missing-path open errors cleanly (no crash).
openErrs:@[{.data.hdb.open["scratch/hdb_eqv_missing_xyz"];0b};(::);{[e] 1b}];

.testutil.assertTrue[pCurveAtTst~hCurveAtTst;"HDB curveAt == parser curveAt (TST)"];
.testutil.assertTrue[pHistTst~hHistTst;"HDB curveHistory == parser curveHistory (TST)"];
.testutil.assertTrue[pCurveAtTst2~hCurveAtTst2;"HDB curveAt == parser curveAt (TST2 - commodity filter isolates)"];
.testutil.assertTrue[hDatesTst~asc distinct pltTst`date;"HDB dates == parser distinct dates (TST)"];
.testutil.assertTrue[openErrs;"open on a missing path errors cleanly"];
.testutil.assertTrue[23=ingRes`rows;"ingested expected row count (TST 4+5+5 + TST2 4+5 = 23)"];

-1 "PASS test_hdb_ingest_query_equivalence: rows=",(string ingRes`rows),", curveAt/curveHistory byte-identical to parser";
