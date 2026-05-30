\l core/init.q
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
bsPrice:.validation.blackScholesClosedForm[`call;100f;100f;1f;0.05;0f;0.2];
optionTable:enlist `optionId`underlying`optionType`strike`expiry`marketPrice!(1;`AAPL;`call;100f;1f;bsPrice);
paramGrid:`initialVarianceList`longRunVarianceList`meanReversionList`volOfVolList`correlationList`jumpIntensityList`jumpMeanList`jumpVolatilityList!(
    enlist 0.04;enlist 0.04;enlist 2.0;enlist 0.0;enlist 0.0;0.0 0.3;enlist 0.0;enlist 0.0);
mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(3000;25;42;0b;0b;0.95);
calibTable:.bates.calibrateGrid[optionTable;mkt;paramGrid;enlist[`mcConfig]!enlist mcConfig];
.testutil.assertTrue[(count calibTable)>0;"calibration has rows"];
bestRow:.bates.bestCalibration calibTable;
.testutil.assertTrue[not null bestRow`rmse;"best rmse finite"];
.testutil.assertTrue[bestRow[`pricedRows]>0;"best has priced rows"];
-1 "PASS test_bates_grid_calibration: gridRows=",string[count calibTable],", bestRmse=",string bestRow`rmse;
