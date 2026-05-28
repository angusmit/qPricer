/ run_benchmarks.q - benchmark runner (separate from unit tests)
\l lib/init.q

-1 "=============================================================================";
-1 " qFDM v",.qfdm.version," Benchmark Suite";
-1 "=============================================================================\n";

benchmarkFiles:enlist "benchmarks/test_benchmark.q";

passCount:0;
failCount:0;

{[benchmarkFile]
    -1 "--- Running: ",benchmarkFile," ---";
    codeLines:read0 hsym `$benchmarkFile;
    filteredLines:codeLines where not codeLines like "\\l *";
    codeBlock:"\n" sv filteredLines;
    testResult:@[{value x;`OK};codeBlock;{x}];
    if[testResult~`OK; passCount+::1];
    if[not testResult~`OK; -2 "  FAILED: ",$[10h=type testResult;testResult;string testResult]; failCount+::1];
    -1 "";
 } each benchmarkFiles;

-1 "=============================================================================";
-1 " Benchmark Results: ",string[passCount]," passed, ",string[failCount]," failed";
-1 "=============================================================================";

if[failCount>0; '"Some benchmarks failed."];
-1 "All benchmarks passed.";
