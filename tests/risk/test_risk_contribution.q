\l core/init.q
/ Build PnL table with group column
pnlTable:();
pnlTable:pnlTable,enlist `scenarioName`bookName`pnl!(`s1;`equities;-5f);
pnlTable:pnlTable,enlist `scenarioName`bookName`pnl!(`s1;`fixed;-3f);
pnlTable:pnlTable,enlist `scenarioName`bookName`pnl!(`s2;`equities;2f);
pnlTable:pnlTable,enlist `scenarioName`bookName`pnl!(`s2;`fixed;1f);
pnlTable:pnlTable,enlist `scenarioName`bookName`pnl!(`s3;`equities;-8f);
pnlTable:pnlTable,enlist `scenarioName`bookName`pnl!(`s3;`fixed;4f);

contrib:.var.riskContribution[pnlTable;`bookName;0.95];
.testutil.assertTrue[2=count contrib;"2 groups"];
firstGroup:contrib 0;
.testutil.assertTrue[`groupVaR in key firstGroup;"has groupVaR"];
.testutil.assertTrue[`totalVaR in key firstGroup;"has totalVaR"];
.testutil.assertTrue[`varContributionPct in key firstGroup;"has contribution"];
.testutil.assertTrue[firstGroup[`status]~`OK;"status OK"];

-1 "PASS test_risk_contribution: groups=",string count contrib;
