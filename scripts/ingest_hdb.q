\l core/init.q
/ ============================================================================
/ scripts/ingest_hdb.q - build the partitioned HDB from the real Barchart CSVs.
/ ----------------------------------------------------------------------------
/ Migration step 3 (ARCHITECTURE.md sections 3 / 9). Reads every commodity folder
/ present under .cfg.paths.barchartRoot via the existing parser (.parser.futures,
/ reused verbatim - expiry/firstDate derived exactly as the parser does) and
/ writes a date-partitioned, splayed HDB at .cfg.paths.hdb. Idempotent: the target
/ is wiped and rebuilt cleanly each run.
/ -
/ This touches REAL, user-supplied, gitignored data and is NOT a test. The suite
/ proves ingest<->query equivalence on SYNTHETIC data instead
/ (tests/data/test_hdb_ingest_query_equivalence.q).
/ -
/ Usage:  q scripts/ingest_hdb.q
/ ============================================================================

barchartRoot:.cfg.paths`barchartRoot;
hdbPath:.cfg.paths`hdb;

/ Discover which commodity folders actually have Day_<TAG>_*.csv files present.
candidates:`CRUDE`GAS`RB`HO;
present:candidates where {[root;tag]
    dir:root,"/",string tag;
    0<count @[{.parser.futures.listFiles[x;y]}[;tag];dir;{[e] ()}]
    }[barchartRoot;] each candidates;

if[0=count present;
    -1 "ingest_hdb: no commodity CSVs found under ",barchartRoot," - nothing to ingest.";
    exit 0];

-1 "ingest_hdb: building HDB at ",hdbPath," from commodities ",(" " sv string present)," ...";
t0:.z.p;
res:.data.hdb.ingest[barchartRoot;present;hdbPath];
elapsed:`long$(.z.p-t0)%1000000;
-1 "ingest_hdb: DONE - ",(string res`rows)," rows across ",(string res`dates),
   " trade dates (splayed) in ",(string elapsed)," ms.";
exit 0;
