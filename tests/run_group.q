/ run_group.q - run a single named test group with timing
/   q tests/run_group.q <groupName> -q
/ Reuses .test.suites from run_all_tests.q without firing its autorun.
.test.skipAutoRun:1b;
\l tests/run_all_tests.q

.test.pass:0;
.test.fail:0;
.test.timings:();
.test.currentGroup:`unknown;

/ Parse group argument
.test.__argv:.z.x;
if[0=count .test.__argv; '"Usage: q tests/run_group.q <groupName> -q"];
.test.__targetGroup:`$first .test.__argv;
.test.__suiteNames:{x 0} each .test.suites;
if[not .test.__targetGroup in .test.__suiteNames;
    '"Unknown group: ",string[.test.__targetGroup]," (available: ",", " sv string .test.__suiteNames,")"];
.test.__matchIdx:.test.__suiteNames?.test.__targetGroup;
.test.__targetSuite:.test.suites .test.__matchIdx;

-1 "=============================================================================";
-1 " qFDM v",.qfdm.version," Group Runner: ",string .test.__targetGroup;
-1 "=============================================================================";

.test.runSuite .test.__targetSuite;

-1 "";
-1 "=============================================================================";
-1 " Group results: ",string[.test.pass]," passed, ",string[.test.fail]," failed";
-1 "=============================================================================";

.test.printTimingReport[];

if[.test.fail>0; '"Group ",string[.test.__targetGroup]," failed: ",string[.test.fail]," failures"];
if[.test.fail=0; -1 "Group passed."];
