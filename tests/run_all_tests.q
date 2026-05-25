/ run_all_tests.q — test runner
/ Usage: q tests/run_all_tests.q

\l lib/init.q

.test.files:(
    "tests/test_european_call.q";
    "tests/test_european_put.q";
    "tests/test_put_call_parity.q";
    "tests/test_grid_convergence.q";
    "tests/test_input_validation.q"
    );

.test.pass:0;
.test.fail:0;

/ Read a test file, strip \l lines, evaluate as one block
.test.run:{[path]
    -1 "--- Running: ",path," ---";
    code:"\n" sv (read0 hsym `$path) where not (read0 hsym `$path) like "\\l *";
    res:@[{value x;`OK};code;{x}];
    if[res~`OK; .test.pass+:1];
    if[not res~`OK; -2 "  FAILED: ",$[10h=type res;res;string res]; .test.fail+:1];
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
