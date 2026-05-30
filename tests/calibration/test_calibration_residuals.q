/ test_calibration_residuals.q
\l core/init.q
spotTable:([] underlying:enlist `AAPL; spot:enlist 100f);
volTable:([] underlying:enlist `AAPL; volatility:enlist 0.2);
rateTable:([] expiry:enlist 1f; riskFreeRate:enlist 0.05);
divTable:([] underlying:enlist `AAPL; dividendYield:enlist 0f);
mktBook:.marketbook.createMarketDataBook[spotTable;volTable;rateTable;divTable];
p1:.validation.blackScholesClosedForm[`call;100f;100f;1f;0.05;0f;0.2];
optionTable:enlist `optionId`underlying`optionType`strike`expiry`marketPrice!(1;`AAPL;`call;100f;1f;p1);

bsCalib:.calibration.calibrateBlackScholesVolSurface[optionTable;mktBook;()!()];
compTable:.modelcompare.priceWithModel[`blackScholesSurface;optionTable;mktBook;bsCalib;()!()];

residualReport:.calibreport.optionResidualReport compTable;
.testutil.assertTrue[1=count residualReport;"1 residual row"];
residualRow:residualReport 0;
.testutil.assertTrue[`residual in key residualRow;"has residual column"];
.testutil.assertTrue[(abs residualRow`residual)<0.001;"residual near zero for BS"];
.testutil.assertNear[residualRow`absoluteError;abs residualRow`residual;0.001;"absError = abs residual"];

-1 "PASS test_calibration_residuals: residual=",string residualRow`residual;
