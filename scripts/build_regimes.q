\l core/init.q
/ ============================================================================
/ scripts/build_regimes.q - build the regimes splay alongside futures in the HDB.
/ ----------------------------------------------------------------------------
/ Research OS R1 (ARCHITECTURE.md 11.1). Opens the real splayed HDB, computes
/ .regime.series for every commodity/date, and writes the `regimes` table splay
/ under .cfg.paths.hdb (mirrors scripts/ingest_hdb.q: .Q.en + p# on commodity + set).
/ Idempotent. Touches real, gitignored data and is NOT a test (the suite proves the
/ regime logic on synthetic data instead). Run AFTER scripts/ingest_hdb.q.
/ -
/ Usage:  q scripts/build_regimes.q
/ ============================================================================
hdbPath:.cfg.paths`hdb;
if[0=count key hsym `$hdbPath,"/sym";
    -1 "build_regimes: no HDB at ",hdbPath," (run scripts/ingest_hdb.q first).";
    exit 0];
.data.hdb.open hdbPath;
commodities:distinct exec commodity from futures;
-1 "build_regimes: labelling regimes for ",(" " sv string commodities)," ...";
t0:.z.p;
res:.regime.buildTable[commodities;hdbPath];
elapsed:`long$(.z.p-t0)%1000000;
-1 "build_regimes: DONE - ",(string res`rows)," regime rows -> ",hdbPath,"/regimes/ in ",(string elapsed)," ms.";
exit 0;
