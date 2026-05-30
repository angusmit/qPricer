/ test_testutil.q - test assertion helpers
\l core/init.q

/ 1. assertTrue passes
.testutil.assertTrue[1b;"true is true"];
.testutil.assertTrue[2>1;"2 > 1"];

/ 2. assertTrue fails correctly
.testutil.expectError["assertTrue false";{.testutil.assertTrue[0b;"should fail"]}];

/ 3. assertEqual passes
.testutil.assertEqual[42;42;"int match"];
.testutil.assertEqual["hello";"hello";"string match"];
.testutil.assertEqual[`sym;`sym;"symbol match"];

/ 4. assertEqual fails correctly
.testutil.expectError["assertEqual mismatch";{.testutil.assertEqual[1;2;"should fail"]}];

/ 5. assertNear passes
.testutil.assertNear[3.14159;3.14;0.01;"pi approx"];

/ 6. assertNear fails correctly
.testutil.expectError["assertNear too far";{.testutil.assertNear[3.14;3.0;0.01;"should fail"]}];

/ 7. assertTableColumns passes
testTbl:([] colA:1 2; colB:3 4);
.testutil.assertTableColumns[testTbl;`colA`colB;"columns present"];

/ 8. assertTableColumns fails correctly
.testutil.expectError["missing column";{.testutil.assertTableColumns[testTbl;`colA`colC;"should fail"]}];

-1 "PASS test_testutil";
