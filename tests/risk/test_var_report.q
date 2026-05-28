\l lib/init.q
pnlVec:-10 -5 -3 -1 0 1 2 3 5 8f;

varRep:.var.varReport[pnlVec;0.9 0.95 0.99];
.testutil.assertTrue[3=count varRep;"3 confidence levels"];

firstRow:varRep 0;
.testutil.assertTrue[`confidenceLevel in key firstRow;"has confidenceLevel"];
.testutil.assertTrue[`valueAtRisk in key firstRow;"has valueAtRisk"];
.testutil.assertTrue[`expectedShortfall in key firstRow;"has expectedShortfall"];
.testutil.assertTrue[`tailCount in key firstRow;"has tailCount"];
.testutil.assertTrue[firstRow[`status]~`OK;"status OK"];

/ From scenario result
scenarioResult:();
idx:0;
pnlList:-10 -5 -3 -1 0 1 2 3 5 8f;
scenNames:`s1`s2`s3`s4`s5`s6`s7`s8`s9`s10;
while[idx<10;
    scenarioResult:scenarioResult,enlist `scenarioName`scenarioPnl!(scenNames idx;pnlList idx);
    idx+:1];
scenRep:.var.varReportFromScenarioResult[scenarioResult;0.95 0.99];
.testutil.assertTrue[2=count scenRep;"2 scenario report rows"];

-1 "PASS test_var_report: rows=",string count varRep;
