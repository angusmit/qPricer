/ test_heston_grid_calibration.q
\l lib/init.q

/ Synthetic BS prices
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
p1:.validation.blackScholesClosedForm[`call;100f;100f;1f;0.05;0f;0.2];
p2:.validation.blackScholesClosedForm[`put;100f;100f;1f;0.05;0f;0.2];
optionTable:();
optionTable:optionTable,enlist `optionId`underlying`optionType`strike`expiry`marketPrice!(1;`AAPL;`call;100f;1f;p1);
optionTable:optionTable,enlist `optionId`underlying`optionType`strike`expiry`marketPrice!(2;`AAPL;`put;100f;1f;p2);

/ Small grid — include xi=0 so BS limit is achievable
paramGridDict:`initialVarianceList`longRunVarianceList`meanReversionList`volOfVolList`correlationList!(
    0.04 0.06;
    0.04 0.06;
    enlist 2.0;
    0.0 0.3;
    -0.7 0.0);

/ Use very small MC for speed
mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(5000;25;42;0b;0b;0.95);
configDict:enlist[`mcConfig]!enlist mcConfig;

gridResult:.calibration.calibrateHestonGrid[optionTable;mkt;paramGridDict;configDict];
.testutil.assertTrue[(count gridResult)>0;"grid has rows"];

bestRow:.calibration.bestCalibrationResult gridResult;
.testutil.assertTrue[not null bestRow`rmse;"best rmse finite"];
.testutil.assertTrue[bestRow[`pricedRows]>0;"best has priced rows"];

/ Grid with xi=0, v0=0.04, theta=0.04 should be among best
-1 "PASS test_heston_grid_calibration: gridRows=",string[count gridResult],", bestRmse=",string[bestRow`rmse],", bestVoV=",string bestRow`volOfVol;
