\l core/init.q

-1 "=============================================================================";
-1 " qFDM v0.31.1 Barchart Historical Option Backtest + Model Comparison";
-1 "=============================================================================";
-1 "";

/ Load and normalise
optionRawTable:.parser.barchart.loadAll "C:/q/w64/qPricer/data/barchart/aapl/options_history";
optionTable:.parser.barchart.normalise optionRawTable;

dataSummary:.parser.barchart.summary optionTable;
-1 "--- Data Summary ---";
-1 "  rows: ",string[dataSummary`rowCount],"  contracts: ",string[dataSummary`contractCount],"  OK: ",string[dataSummary`okRows],"  errors: ",string dataSummary`errorRows;
-1 "  dates: ",string[dataSummary`minDate]," to ",string dataSummary`maxDate;
-1 "";

/ Model pricing
pricingConfig:`riskFreeRate`dividendYield`pricingModel!(0.05;0.005;`blackScholes);
pricedOptionTable:.parser.barchart.priceOptionRows[optionTable;pricingConfig];
pricedOkCnt:sum (pricedOptionTable`pricingStatus)=`OK;
-1 "--- Model Pricing ---";
-1 "  priced OK: ",string[pricedOkCnt]," of ",string count pricedOptionTable;
-1 "";

/ Multi-day replay
availDates:.parser.barchart.availableDates optionTable;
startDt:first availDates;
endDt:last availDates;
replay:.parser.barchart.multiDayReplay[optionTable;startDt;endDt;1f];

bSummary:.parser.barchart.backtestSummary replay;
-1 "--- Backtest Summary ---";
-1 "  rows: ",string[bSummary`rowCount],"  contracts: ",string[bSummary`contractCount],"  errors: ",string bSummary`errorRows;
-1 "  period: ",string[bSummary`startDate]," to ",string bSummary`endDate;
-1 "  totalMarketPnl: ",string bSummary`totalMarketPnl;
-1 "";

/ Model replay
modelReplayResult:.parser.barchart.modelReplay[replay;pricedOptionTable];
compSummary:.parser.barchart.modelComparisonSummary modelReplayResult;
-1 "--- Model vs Market Summary ---";
-1 "  totalMarketPnl: ",string compSummary`totalMarketPnl;
-1 "  totalModelPnl: ",string compSummary`totalModelPnl;
-1 "  totalPnlDiff: ",string compSummary`totalPnlDifference;
-1 "  meanAbsEntryError: ",string compSummary`meanAbsEntryModelError;
-1 "  meanAbsExitError: ",string compSummary`meanAbsExitModelError;
-1 "  meanAbsPnlDiff: ",string compSummary`meanAbsPnlDifference;
-1 "";

/ Error buckets
mBuckets:.parser.barchart.errorByMoneyness pricedOptionTable;
-1 "--- Pricing Error by Moneyness ---";
mbIdx:0;
while[mbIdx<count mBuckets;
    bRow:mBuckets mbIdx;
    -1 "  ",string[bRow`bucket]," n=",string[bRow`rowCount]," meanErr=",string[bRow`meanModelError]," meanAbsErr=",string bRow`meanAbsModelError;
    mbIdx+:1];
-1 "";

dBuckets:.parser.barchart.errorByDte pricedOptionTable;
-1 "--- Pricing Error by DTE ---";
dbIdx:0;
while[dbIdx<count dBuckets;
    bRow:dBuckets dbIdx;
    -1 "  ",string[bRow`bucket]," n=",string[bRow`rowCount]," meanErr=",string[bRow`meanModelError]," meanAbsErr=",string bRow`meanAbsModelError;
    dbIdx+:1];
-1 "";

oBuckets:.parser.barchart.errorByOptionType pricedOptionTable;
-1 "--- Pricing Error by Option Type ---";
obIdx:0;
while[obIdx<count oBuckets;
    bRow:oBuckets obIdx;
    -1 "  ",string[bRow`bucket]," n=",string[bRow`rowCount]," meanErr=",string[bRow`meanModelError]," meanAbsErr=",string bRow`meanAbsModelError;
    obIdx+:1];
-1 "";

/ Performance
daily:.parser.barchart.dailyPnl replay;
perf:.parser.barchart.performanceSummary daily;
-1 "--- Performance Summary ---";
-1 "  observations: ",string perf`observationCount;
-1 "  totalPnl: ",string perf`totalMarketPnl;
-1 "  meanDaily: ",string perf`meanDailyPnl;
-1 "  dailyVol: ",string perf`dailyPnlVolatility;
-1 "  winRate: ",string perf`winRate;
-1 "  maxDrawdown: ",string perf`maxDrawdown;
-1 "";

-1 "=============================================================================";
-1 " Barchart backtest + model comparison completed";
-1 "=============================================================================";
