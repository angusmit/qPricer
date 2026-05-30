\l core/init.q
/ Build scenario result from synthetic data
scenarioResult:();
scenNames:`spotUp5`spotUp2`base`spotDn2`spotDn5`volUp`volDn`rateUp`rateDn;
pnlValues:6.5 2.5 0.1 -2.3 -6.1 1.2 -0.8 0.5 -0.3;
sIdx:0;
while[sIdx<count scenNames;
    scenarioResult:scenarioResult,enlist `scenarioName`scenarioPnl!(scenNames sIdx;pnlValues sIdx);
    sIdx+:1];

/ Build PnL distribution
pnlDist:.riskdist.buildPnlDistribution scenarioResult;
pnlVec:.riskdist.pnlVector pnlDist;
.testutil.assertTrue[9=count pnlVec;"9 scenario PnLs"];

/ VaR report
varRep:.var.varReport[pnlVec;0.95 0.99];
.testutil.assertTrue[2=count varRep;"2 confidence levels"];
row95:varRep 0;
.testutil.assertTrue[not null row95`valueAtRisk;"VaR finite"];
.testutil.assertTrue[not null row95`expectedShortfall;"ES finite"];
.testutil.assertTrue[row95[`status]~`OK;"status OK"];

-1 "PASS test_portfolio_var: VaR95=",string[row95`valueAtRisk],", ES95=",string row95`expectedShortfall;
