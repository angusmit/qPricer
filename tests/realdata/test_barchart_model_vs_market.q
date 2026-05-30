\l core/init.q
/ Build synthetic priced option table and replay for model comparison
pricingConfig:`riskFreeRate`dividendYield`pricingModel!(0.05;0.005;`blackScholes);

/ 2 dates, 1 contract
optionTable:();
optionTable:optionTable,enlist `snapshotDate`underlying`contractId`expiryDate`optionType`strike`spot`open`high`low`latest`bid`ask`mid`marketPrice`impliedVolatility`volume`openInterest`vendorDelta`vendorGamma`vendorTheta`vendorVega`vendorRho`vendorTheo`dte`tau`moneyness`bidAskSpread`bidAskSpreadPct`status`errorMessage!(
    2024.01.02;`aapl;`testC1;2024.02.16;`call;165f;185f;22f;23f;21f;22.05;21.90;22.20;22.05;22.05;0.2755;100f;2000f;0.91;0.009;-0.053;0.104;0.177;22.4;45i;45%365f;165%185f;0.30;0.014;`OK;"");
optionTable:optionTable,enlist `snapshotDate`underlying`contractId`expiryDate`optionType`strike`spot`open`high`low`latest`bid`ask`mid`marketPrice`impliedVolatility`volume`openInterest`vendorDelta`vendorGamma`vendorTheta`vendorVega`vendorRho`vendorTheo`dte`tau`moneyness`bidAskSpread`bidAskSpreadPct`status`errorMessage!(
    2024.01.03;`aapl;`testC1;2024.02.16;`call;165f;186f;23f;24f;22f;23.10;22.95;23.25;23.10;23.10;0.28;100f;2000f;0.92;0.008;-0.052;0.103;0.178;23.5;44i;44%365f;165%186f;0.30;0.013;`OK;"");

pricedTable:.parser.barchart.priceOptionRows[optionTable;pricingConfig];
.testutil.assertTrue[2=count pricedTable;"2 priced rows"];

/ Build replay
replay:.parser.barchart.oneDayReplay[optionTable;2024.01.02;2024.01.03;1f];
.testutil.assertTrue[1=count replay;"1 replay row"];

/ Model replay
modelReplayResult:.parser.barchart.modelReplay[replay;pricedTable];
.testutil.assertTrue[1=count modelReplayResult;"1 model replay row"];
mRow:modelReplayResult 0;
.testutil.assertTrue[not null mRow`entryModelPrice;"entry model price populated"];
.testutil.assertTrue[not null mRow`exitModelPrice;"exit model price populated"];

/ modelPnl = qty * 100 * (exit - entry)
expectedModelPnl:1f*100f*mRow[`exitModelPrice]-mRow`entryModelPrice;
.testutil.assertNear[mRow`modelPnl;expectedModelPnl;0.01;"modelPnl formula"];

/ pnlDifference = modelPnl - marketPnl
.testutil.assertNear[mRow`pnlDifference;mRow[`modelPnl]-mRow`marketPnl;0.01;"pnlDiff formula"];

/ Comparison summary
compSummary:.parser.barchart.modelComparisonSummary modelReplayResult;
.testutil.assertTrue[compSummary[`rowCount]=1;"1 comparison row"];
.testutil.assertTrue[not null compSummary`totalModelPnl;"totalModelPnl populated"];

/ Error buckets (at least return rows without crash)
mBuckets:.parser.barchart.errorByMoneyness pricedTable;
.testutil.assertTrue[(count mBuckets)>0;"moneyness buckets exist"];
dBuckets:.parser.barchart.errorByDte pricedTable;
.testutil.assertTrue[(count dBuckets)>0;"dte buckets exist"];
oBuckets:.parser.barchart.errorByOptionType pricedTable;
.testutil.assertTrue[(count oBuckets)>0;"option type buckets exist"];

-1 "PASS test_barchart_model_vs_market: modelPnl=",string[mRow`modelPnl],", pnlDiff=",string mRow`pnlDifference;
