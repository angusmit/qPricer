/ test_merton_model_comparison.q
\l core/init.q
spotTable:([] underlying:enlist `AAPL; spot:enlist 100f);
volTable:([] underlying:enlist `AAPL; volatility:enlist 0.2);
rateTable:([] expiry:enlist 1f; riskFreeRate:enlist 0.05);
divTable:([] underlying:enlist `AAPL; dividendYield:enlist 0f);
mktBook:.marketbook.createMarketDataBook[spotTable;volTable;rateTable;divTable];

bsPrice:.validation.blackScholesClosedForm[`call;100f;100f;1f;0.05;0f;0.2];
optionTable:enlist `optionId`underlying`optionType`strike`expiry`marketPrice!(1;`AAPL;`call;100f;1f;bsPrice);

/ BS calibration
bsCalib:.calibration.calibrateBlackScholesVolSurface[optionTable;mktBook;()!()];

/ Merton calibration result (lambda=0 should match BS)
mertonCalib:`volatility`jumpIntensity`jumpMean`jumpVolatility!(0.2;0.0;0.0;0.0);

calibDict:`blackScholesSurface`mertonSeries!(bsCalib;mertonCalib);
compTable:.modelcompare.compareModels[`blackScholesSurface`mertonSeries;optionTable;mktBook;calibDict;()!()];

.testutil.assertTrue[(count compTable)>0;"comparison has rows"];
ranked:.modelcompare.rankModels compTable;
.testutil.assertTrue[2=count ranked;"2 models ranked"];

-1 "PASS test_merton_model_comparison: rows=",string count compTable;
