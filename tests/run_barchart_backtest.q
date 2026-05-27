\l lib/init.q

-1 "=============================================================================";
-1 " qFDM v0.31 Barchart Historical Option Backtest";
-1 "=============================================================================";
-1 "";

/ Load raw data
optionRawTable:.parser.barchart.loadAll "C:/q/w64/qPricer/data/barchart/aapl/options_history";
optionTable:.parser.barchart.normalise optionRawTable;

/ Data summary
dataSummary:.parser.barchart.summary optionTable;
-1 "--- Data Summary ---";
-1 "  rows: ",string[dataSummary`rowCount],"  contracts: ",string[dataSummary`contractCount],"  OK: ",string[dataSummary`okRows],"  errors: ",string dataSummary`errorRows;
-1 "  dates: ",string[dataSummary`minDate]," to ",string dataSummary`maxDate;
-1 "";

/ Multi-day replay
availDates:.parser.barchart.availableDates optionTable;
startDt:first availDates;
endDt:last availDates;
-1 "Running replay from ",string[startDt]," to ",string[endDt],"...";
replay:.parser.barchart.multiDayReplay[optionTable;startDt;endDt;1f];

/ Backtest summary
bSummary:.parser.barchart.backtestSummary replay;
-1 "";
-1 "--- Backtest Summary ---";
-1 "  rows: ",string[bSummary`rowCount],"  contracts: ",string[bSummary`contractCount],"  errors: ",string bSummary`errorRows;
-1 "  period: ",string[bSummary`startDate]," to ",string bSummary`endDate;
-1 "  totalPnl: ",string[bSummary`totalMarketPnl],"  worst: ",string[bSummary`worstTradePnl],"  best: ",string bSummary`bestTradePnl;
-1 "";

/ Daily PnL
daily:.parser.barchart.dailyPnl replay;
-1 "--- Daily PnL (first 10) ---";
showIdx:0;
showCnt:10&count daily;
while[showIdx<showCnt;
    dRow:daily showIdx;
    -1 "  ",string[dRow`snapshotDate],"  daily: ",string[dRow`dailyMarketPnl],"  cum: ",string[dRow`cumulativeMarketPnl],"  trades: ",string dRow`tradeCount;
    showIdx+:1];
if[showCnt<count daily; -1 "  ... (",string[(count daily)-showCnt]," more rows)"];
-1 "";

/ Performance summary
perf:.parser.barchart.performanceSummary daily;
-1 "--- Performance Summary ---";
-1 "  observations: ",string perf`observationCount;
-1 "  totalPnl: ",string perf`totalMarketPnl;
-1 "  meanDaily: ",string perf`meanDailyPnl;
-1 "  dailyVol: ",string perf`dailyPnlVolatility;
-1 "  winRate: ",string perf`winRate;
-1 "  worstDay: ",string perf`worstDailyPnl;
-1 "  bestDay: ",string perf`bestDailyPnl;
-1 "  maxDrawdown: ",string perf`maxDrawdown;
-1 "";

-1 "=============================================================================";
-1 " Barchart backtest completed";
-1 "=============================================================================";
