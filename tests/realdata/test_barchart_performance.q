\l lib/init.q
/ Build simple dailyPnl table
dailyPnlTable:();
dailyPnlTable:dailyPnlTable,enlist `snapshotDate`dailyMarketPnl`cumulativeMarketPnl`tradeCount`okRows`errorRows`status`errorMessage!(2024.01.03;500f;500f;4;4;0;`OK;"");
dailyPnlTable:dailyPnlTable,enlist `snapshotDate`dailyMarketPnl`cumulativeMarketPnl`tradeCount`okRows`errorRows`status`errorMessage!(2024.01.04;-200f;300f;4;4;0;`OK;"");
dailyPnlTable:dailyPnlTable,enlist `snapshotDate`dailyMarketPnl`cumulativeMarketPnl`tradeCount`okRows`errorRows`status`errorMessage!(2024.01.05;100f;400f;4;4;0;`OK;"");

/ Performance summary
perf:.parser.barchart.performanceSummary dailyPnlTable;
.testutil.assertTrue[perf[`observationCount]=3;"3 observations"];
.testutil.assertNear[perf`totalMarketPnl;400f;0.01;"total PnL = 400"];
.testutil.assertTrue[(perf`winRate)>0f;"winRate > 0"];
.testutil.assertTrue[(perf`winRate)<=1f;"winRate <= 1"];
.testutil.assertTrue[(perf`maxDrawdown)>=0f;"maxDrawdown >= 0"];
.testutil.assertNear[perf`maxDrawdown;200f;0.01;"maxDrawdown = 200 (500 peak to 300)"];

/ Drawdown
ddTable:.parser.barchart.drawdown dailyPnlTable;
.testutil.assertTrue[3=count ddTable;"3 drawdown rows"];
.testutil.assertTrue[`runningPeak in cols ddTable;"has runningPeak"];
.testutil.assertTrue[`drawdown in cols ddTable;"has drawdown"];
.testutil.assertNear[(ddTable 1)[`drawdown];200f;0.01;"day2 drawdown=200"];

-1 "PASS test_barchart_performance: maxDrawdown=",string[perf`maxDrawdown],", winRate=",string perf`winRate;
