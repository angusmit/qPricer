/ test_model_comparison_pricing.q
\l core/init.q
spotTable:([] underlying:enlist `AAPL; spot:enlist 100f);
volTable:([] underlying:enlist `AAPL; volatility:enlist 0.2);
rateTable:([] expiry:enlist 1f; riskFreeRate:enlist 0.05);
divTable:([] underlying:enlist `AAPL; dividendYield:enlist 0f);
mktBook:.marketbook.createMarketDataBook[spotTable;volTable;rateTable;divTable];
bsPrice:.validation.blackScholesClosedForm[`call;100f;100f;1f;0.05;0f;0.2];
optionTable:enlist `optionId`underlying`optionType`strike`expiry`marketPrice!(1;`AAPL;`call;100f;1f;bsPrice);

bsCalib:.calibration.calibrateBlackScholesVolSurface[optionTable;mktBook;()!()];
compTable:.modelcompare.priceWithModel[`blackScholesSurface;optionTable;mktBook;bsCalib;()!()];

calibDict:enlist[`blackScholesSurface]!enlist bsCalib;
tradeTable:enlist `tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;100000f);

pricingResult:.modelcompare.applyBestModel[tradeTable;mktBook;compTable;calibDict;()!()];
resultRow:pricingResult 0;
.testutil.assertTrue[resultRow[`status]~`OK;"pricing OK"];
.testutil.assertTrue[resultRow[`unitPrice]>0f;"price positive"];
.testutil.assertNear[resultRow`unitPrice;bsPrice;0.01;"price matches BS"];

-1 "PASS test_model_comparison_pricing: price=",string resultRow`unitPrice;
