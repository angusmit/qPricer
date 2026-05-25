/ run_smoke_tests.q - fast core checks (~10 seconds)
\l ./lib/init.q

.test.files:(
    "tests/test_european_call.q";
    "tests/test_european_put.q";
    "tests/test_input_validation.q";
    "tests/test_implied_vol_call.q";
    "tests/test_result.q";
    "tests/test_testutil.q"
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
-1 " qFDM v",.qfdm.version," Smoke Tests";
-1 "=============================================================================\n";

.test.run each .test.files;

-1 "=============================================================================";
-1 " Smoke Results: ",string[.test.pass]," passed, ",string[.test.fail]," failed";
-1 "=============================================================================";

if[.test.fail=0; -1 "All smoke tests passed."];
if[.test.fail>0; '"Some smoke tests failed."];
