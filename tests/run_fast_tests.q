/ run_fast_tests.q - fast-tier runner (deterministic, no heavy MC, coarse grids)
/ Target: < 2 minutes wall time. Used for quick verification during development.
.test.skipAutoRun:1b;
\l tests/run_all_tests.q

.test.pass:0;
.test.fail:0;
.test.timings:();
.test.currentGroup:`unknown;

/ Fast tier composition: core + greeks + market + impliedvol + infra + a subset of
/ deterministic strategy/commodity tests. Excludes heavy Monte Carlo and large FDM
/ runs.
.test.fastSuites:(
    (`core;       .test.coreFiles);
    (`greeks;     .test.greeksFiles);
    (`market;     .test.marketFiles);
    (`impliedvol; .test.impliedVolFiles);
    (`infra;      .test.infraFiles));

-1 "=============================================================================";
-1 " qFDM v",.qfdm.version," Fast Tier";
-1 "=============================================================================";

.test.runSuite each .test.fastSuites;

-1 "";
-1 "=============================================================================";
-1 " Fast results: ",string[.test.pass]," passed, ",string[.test.fail]," failed";
-1 "=============================================================================";

.test.printTimingReport[];

if[.test.fail>0; '"Fast suite failed: ",string[.test.fail]," failures"];
if[.test.fail=0; -1 "Fast suite passed."];
