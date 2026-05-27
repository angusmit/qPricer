\l lib/init.q
spotTable:([] underlying:enlist `AAPL; spot:enlist 100f);
volTable:([] underlying:enlist `AAPL; volatility:enlist 0.2);
rateTable:([] expiry:enlist 1f; riskFreeRate:enlist 0.05);
divTable:([] underlying:enlist `AAPL; dividendYield:enlist 0f);
mktBook:.marketbook.createMarketDataBook[spotTable;volTable;rateTable;divTable];
bsPrice:.validation.blackScholesClosedForm[`call;100f;100f;1f;0.05;0f;0.2];
optionTable:enlist `optionId`underlying`optionType`strike`expiry`marketPrice!(1;`AAPL;`call;100f;1f;bsPrice);
bsCalib:.calibration.calibrateBlackScholesVolSurface[optionTable;mktBook;()!()];
batesCalib:`initialVariance`longRunVariance`meanReversion`volOfVol`correlation`jumpIntensity`jumpMean`jumpVolatility!(0.04;0.04;2.0;0.0;0.0;0.0;0.0;0.0);
mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(3000;25;42;0b;0b;0.95);
calibDict:`blackScholesSurface`bates!(bsCalib;batesCalib);
compTable:.modelcompare.compareModels[`blackScholesSurface`bates;optionTable;mktBook;calibDict;enlist[`mcConfig]!enlist mcConfig];
.testutil.assertTrue[(count compTable)>0;"comparison has rows"];
ranked:.modelcompare.rankModels compTable;
.testutil.assertTrue[2=count ranked;"2 models ranked"];
-1 "PASS test_bates_model_comparison: rows=",string count compTable;
