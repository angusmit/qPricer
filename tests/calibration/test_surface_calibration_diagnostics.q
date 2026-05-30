/ test_surface_calibration_diagnostics.q
\l core/init.q

spotTable:([] underlying:enlist `AAPL; spot:enlist 100f);
volTable:([] underlying:enlist `AAPL; volatility:enlist 0.2);
rateTable:([] expiry:enlist 1f; riskFreeRate:enlist 0.05);
divTable:([] underlying:enlist `AAPL; dividendYield:enlist 0f);
mktBook:.marketbook.createMarketDataBook[spotTable;volTable;rateTable;divTable];
p1:.validation.blackScholesClosedForm[`call;100f;100f;1f;0.05;0f;0.2];
optionTable:enlist `optionId`underlying`optionType`strike`expiry`marketPrice!(1;`AAPL;`call;100f;1f;p1);

calibResult:.calibration.calibrateBlackScholesVolSurface[optionTable;mktBook;()!()];

surfaceTable:.calibration.buildVolSurfaceFromCalibration calibResult;
.testutil.assertTrue[1=count surfaceTable;"1 surface row"];
surfaceRow:surfaceTable 0;
.testutil.assertNear[surfaceRow`impliedVolatility;0.2;0.001;"surface IV=0.2"];

diagResult:.calibration.surfaceDiagnostics calibResult;
.testutil.assertTrue[diagResult[`optionCount]=1;"diag optionCount"];
.testutil.assertNear[diagResult`averageIV;0.2;0.001;"diag avgIV"];
.testutil.assertTrue[diagResult[`failedRows]=0;"diag no failures"];

-1 "PASS test_surface_calibration_diagnostics: avgIV=",string diagResult`averageIV;
