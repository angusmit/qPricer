/ test_calibrated_model_pricing.q
\l core/init.q

/ Calibrate BS vol surface
spotTable:([] underlying:enlist `AAPL; spot:enlist 100f);
volTable:([] underlying:enlist `AAPL; volatility:enlist 0.2);
rateTable:([] expiry:enlist 1f; riskFreeRate:enlist 0.05);
divTable:([] underlying:enlist `AAPL; dividendYield:enlist 0f);
mktBook:.marketbook.createMarketDataBook[spotTable;volTable;rateTable;divTable];

bsPrice:.validation.blackScholesClosedForm[`call;100f;100f;1f;0.05;0f;0.2];
optionTable:enlist `optionId`underlying`optionType`strike`expiry`marketPrice!(1;`AAPL;`call;100f;1f;bsPrice);
calibResult:.calibration.calibrateBlackScholesVolSurface[optionTable;mktBook;()!()];

/ Use the first calibrated row as model input
calibRow:calibResult 0;
/ Price a trade using calibrated IV
tradeTable:enlist `tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;100000f);
pricingResult:.calibration.applyCalibratedModel[tradeTable;mktBook;calibRow;()!()];

resultRow:pricingResult 0;
.testutil.assertTrue[resultRow[`status]~`OK;"calibrated pricing OK"];
.testutil.assertTrue[resultRow[`unitPrice]>0f;"calibrated price positive"];
.testutil.assertNear[resultRow`unitPrice;bsPrice;0.01;"calibrated price matches BS"];
.testutil.assertNear[resultRow`notionalPrice;resultRow[`unitPrice]*100000f;1f;"notional correct"];

-1 "PASS test_calibrated_model_pricing: price=",string[resultRow`unitPrice],", expected=",string bsPrice;
