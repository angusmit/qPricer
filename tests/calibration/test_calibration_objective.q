/ test_calibration_objective.q
\l core/init.q

.testutil.assertNear[.objective.squaredError[10f;9f];1f;0.001;"squaredError"];
.testutil.assertNear[.objective.absoluteError[10f;9f];1f;0.001;"absoluteError"];
.testutil.assertNear[.objective.relativeError[10f;9f];1f%9f;0.001;"relativeError"];

modelVec:10 11 12f;
mktVec:10.1 10.8 12.3f;
.testutil.assertTrue[.objective.rmse[modelVec;mktVec]>0f;"rmse positive"];
.testutil.assertTrue[.objective.mae[modelVec;mktVec]>0f;"mae positive"];

wsse:.objective.weightedSse[modelVec;mktVec;1 1 1f];
.testutil.assertTrue[wsse>0f;"weightedSse positive"];

/ Summary
resultTable:();
resultTable:resultTable,enlist `modelPrice`marketPrice`absoluteError!(10f;10.1;0.1);
resultTable:resultTable,enlist `modelPrice`marketPrice`absoluteError!(11f;10.8;0.2);
summaryResult:.objective.calibrationSummary resultTable;
.testutil.assertTrue[summaryResult[`rowCount]=2;"summary rowCount"];
.testutil.assertTrue[summaryResult[`rmse]>0f;"summary rmse"];
.testutil.assertTrue[summaryResult[`mae]>0f;"summary mae"];

-1 "PASS test_calibration_objective";
