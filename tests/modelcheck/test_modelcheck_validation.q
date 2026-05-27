\l lib/init.q

/ Empty summary works
emptySummary:.limitcheck.summary ();
.testutil.assertTrue[emptySummary[`checkCount]=0;"empty summary checkCount=0"];

/ resultRow works
testRow:.limitcheck.resultRow["test";`a;`b;10f;10.1;0.5;0.1];
.testutil.assertTrue[testRow`passed;"within tolerance"];
.testutil.assertNear[testRow`absoluteDifference;0.1;0.001;"abs diff correct"];

/ errorRow works
errRow:.limitcheck.errorRow["test";`a;`b;"some error"];
.testutil.assertTrue[errRow[`status]~`ERROR;"error status"];
.testutil.assertTrue[not errRow`passed;"error not passed"];

-1 "PASS test_modelcheck_validation";
