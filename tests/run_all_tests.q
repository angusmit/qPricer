/ run_all_tests.q — test runner
\l lib/init.q

.test.files:(
    "tests/test_european_call.q";
    "tests/test_european_put.q";
    "tests/test_put_call_parity.q";
    "tests/test_grid_convergence.q";
    "tests/test_input_validation.q";
    "tests/test_greeks_call.q";
    "tests/test_greeks_put.q";
    "tests/test_scenario_risk.q";
    "tests/test_american_put.q"
    );

.test.pass:0;
.test.fail:0;

.test.run:{[testPath]
    -1 "--- Running: ",testPath," ---";
    codeLines:read0 hsym `$testPath;
    filteredLines:codeLines where not codeLines like "\\l *";
    codeBlock:"\n" sv filteredLines;
    testResult:@[{value x;`OK};codeBlock;{x}];
    if[testResult~`OK; .test.pass+:1];
    if[not testResult~`OK; -2 "  FAILED: ",$[10h=type testResult;testResult;string testResult]; .test.fail+:1];
    -1 "";
 };

-1 "=============================================================================";
-1 " qFDM v",.qfdm.version," Test Suite";
-1 "=============================================================================\n";

.test.run each .test.files;

-1 "=============================================================================";
-1 " Results: ",string[.test.pass]," passed, ",string[.test.fail]," failed";
-1 "=============================================================================";

if[.test.fail=0; -1 "All tests passed."];
if[.test.fail>0; '"Some tests failed: ",string[.test.fail]," failures"];
