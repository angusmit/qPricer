/ test_black_scholes_vol_calibration.q
\l lib/init.q

spotTable:([] underlying:enlist `AAPL; spot:enlist 100f);
volTable:([] underlying:enlist `AAPL; volatility:enlist 0.2);
rateTable:([] expiry:enlist 1f; riskFreeRate:enlist 0.05);
divTable:([] underlying:enlist `AAPL; dividendYield:enlist 0f);
mktBook:.marketbook.createMarketDataBook[spotTable;volTable;rateTable;divTable];

p1:.validation.blackScholesClosedForm[`call;100f;90f;0.5;0.05;0f;0.2];
p2:.validation.blackScholesClosedForm[`call;100f;100f;1f;0.05;0f;0.2];
p3:.validation.blackScholesClosedForm[`put;100f;110f;1f;0.05;0f;0.2];

optionTable:();
optionTable:optionTable,enlist `optionId`underlying`optionType`strike`expiry`marketPrice!(1;`AAPL;`call;90f;0.5;p1);
optionTable:optionTable,enlist `optionId`underlying`optionType`strike`expiry`marketPrice!(2;`AAPL;`call;100f;1f;p2);
optionTable:optionTable,enlist `optionId`underlying`optionType`strike`expiry`marketPrice!(3;`AAPL;`put;110f;1f;p3);

calibResult:.calibration.calibrateBlackScholesVolSurface[optionTable;mktBook;()!()];

/ Use column access (calibResult is a table)
okCount:sum calibResult[`status]=`OK;
.testutil.assertTrue[okCount=3;"all 3 rows calibrated"];

ivs:calibResult`impliedVolatility;
avgIV:avg ivs;
.testutil.assertNear[avgIV;0.2;0.001;"average IV close to 0.2"];

maxErr:max calibResult`absoluteError;
.testutil.assertTrue[maxErr<0.001;"model prices match"];

-1 "PASS test_black_scholes_vol_calibration: avgIV=",string[avgIV],", maxErr=",string maxErr;
