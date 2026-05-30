/ test_sabr_model_comparison.q
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

/ SABR calibration (best row serves as calibration result)
sabrBest:`alpha`beta`rho`nu!(0.2;0.5;0.0;0.0001);

/ Compare models
calibDict:`blackScholesSurface`sabrSurface!(bsCalib;sabrBest);
compTable:.modelcompare.compareModels[`blackScholesSurface`sabrSurface;optionTable;mktBook;calibDict;()!()];

.testutil.assertTrue[(count compTable)>0;"comparison has rows"];
ranked:.modelcompare.rankModels compTable;
.testutil.assertTrue[2=count ranked;"2 models ranked"];
bestRow:.modelcompare.bestModel compTable;
.testutil.assertTrue[not null bestRow`rmse;"best rmse finite"];

-1 "PASS test_sabr_model_comparison: best=",string[bestRow`modelName],", rmse=",string bestRow`rmse;
