/ test_result.q - standardised result utilities
\l core/init.q

/ 1. okRow
okDict:.result.okRow[`tradeId`unitPrice!(1;10.5)];
.testutil.assertEqual[okDict`status;`OK;"okRow status"];
.testutil.assertEqual[okDict`errorMessage;"";"okRow errorMessage"];
.testutil.assertEqual[okDict`unitPrice;10.5;"okRow preserves data"];

/ 2. errorRow
errDict:.result.errorRow[`tradeId`unitPrice!(2;0Nf);"price failed"];
.testutil.assertEqual[errDict`status;`ERROR;"errorRow status"];
.testutil.assertTrue[0<count errDict`errorMessage;"errorRow has message"];

/ 3. isOk / isError on table
resultTbl:();
resultTbl:resultTbl,enlist .result.okRow[`tradeId`unitPrice!(1;10.5)];
resultTbl:resultTbl,enlist .result.errorRow[`tradeId`unitPrice!(2;0Nf);"failed"];
okMask:.result.isOk resultTbl;
.testutil.assertEqual[okMask;10b;"isOk mask"];
errorMask:.result.isError resultTbl;
.testutil.assertEqual[errorMask;01b;"isError mask"];

/ 4. errorSummary
errorRows:.result.errorSummary resultTbl;
.testutil.assertEqual[count errorRows;1;"errorSummary count"];

/ 5. assertColumns
.result.assertColumns[resultTbl;`tradeId`unitPrice`status`errorMessage];
.testutil.expectError["assertColumns missing col";{.result.assertColumns[([]a:1 2);`a`b`c]}];

-1 "PASS test_result";
