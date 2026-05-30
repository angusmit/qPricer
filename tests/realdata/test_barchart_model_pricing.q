\l core/init.q
/ Synthetic option row for BS pricing test
/ AAPL call: S=185, K=165, tau=45/365, vol=0.2755, r=0.05, q=0.005
pricingConfig:`riskFreeRate`dividendYield`pricingModel!(0.05;0.005;`blackScholes);

optionTable:();
optionTable:optionTable,enlist `snapshotDate`underlying`contractId`expiryDate`optionType`strike`spot`open`high`low`latest`bid`ask`mid`marketPrice`impliedVolatility`volume`openInterest`vendorDelta`vendorGamma`vendorTheta`vendorVega`vendorRho`vendorTheo`dte`tau`moneyness`bidAskSpread`bidAskSpreadPct`status`errorMessage!(
    2024.01.02;`aapl;`testCall;2024.02.16;`call;165f;185f;22f;23f;21f;22.05;21.90;22.20;22.05;22.05;0.2755;100f;2000f;0.91;0.009;-0.053;0.104;0.177;22.4;45i;45%365f;165f%185f;0.30;0.014;`OK;"");
optionTable:optionTable,enlist `snapshotDate`underlying`contractId`expiryDate`optionType`strike`spot`open`high`low`latest`bid`ask`mid`marketPrice`impliedVolatility`volume`openInterest`vendorDelta`vendorGamma`vendorTheta`vendorVega`vendorRho`vendorTheo`dte`tau`moneyness`bidAskSpread`bidAskSpreadPct`status`errorMessage!(
    2024.01.02;`aapl;`testPut;2024.02.16;`put;185f;185f;3.50;3.80;3.20;3.50;3.40;3.60;3.50;3.50;0.30;200f;1500f;-0.25;0.015;-0.040;0.080;-0.100;3.6;45i;45%365f;1f;0.20;0.057;`OK;"");

pricedTable:.parser.barchart.priceOptionRows[optionTable;pricingConfig];

/ Call model price should be positive and reasonable (deep ITM ~20)
callRow:pricedTable 0;
.testutil.assertTrue[callRow[`pricingStatus]=`OK;"call pricing OK"];
.testutil.assertTrue[callRow[`modelPrice]>0f;"call modelPrice positive"];
.testutil.assertTrue[callRow[`modelPrice]>15f;"call modelPrice > 15 (deep ITM)"];
.testutil.assertTrue[callRow[`modelPrice]<30f;"call modelPrice < 30"];

/ Put model price should be positive
putRow:pricedTable 1;
.testutil.assertTrue[putRow[`pricingStatus]=`OK;"put pricing OK"];
.testutil.assertTrue[putRow[`modelPrice]>0f;"put modelPrice positive"];

/ Model error = modelPrice - marketPrice
.testutil.assertNear[callRow`modelError;callRow[`modelPrice]-callRow`marketPrice;0.001;"modelError formula"];

/ Greeks populated
.testutil.assertTrue[not null callRow`modelDelta;"call delta populated"];
.testutil.assertTrue[callRow[`modelDelta]>0.5;"call delta > 0.5 (ITM)"];

-1 "PASS test_barchart_model_pricing: callModel=",string[callRow`modelPrice],", putModel=",string putRow`modelPrice;
