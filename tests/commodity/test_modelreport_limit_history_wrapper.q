\l lib/init.q
/ runPortfolioRiskWithLimitHistory wraps a single portfolio run, generates
/ snapshots, appends to existing histories, and rebuilds the dashboard. Verify
/ histories grow across two successive calls. Keep MC paths small for speed.

mcCfg:.montecarlo.defaultMcConfig[];
mcCfg:@[mcCfg;`pathCount;:;2000];
mcCfg:@[mcCfg;`timeStepCount;:;15];

setupCall:`optionType`strikePrice`expiry`riskFreeRate`forwardPrice`spotPrice!(`call;75f;1f;0.05;75f;75f);
inputsAll:`black76`schwartz`schwartz2`mrjump!(
    enlist[`volatility]!enlist 0.30;
    `x0`params!(log 75f;`meanReversionSpeed`longRunLogMean`volatility!(1.5;log 75f;0.30));
    `shortFactor0`longFactor0`params!(0f;log 75f;`meanReversionSpeed`shortVolatility`longVolatility`longDrift`correlation!(1.5;0.30;0.15;0f;0.25));
    `x0`params`mcConfig!(log 75f;`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility!(2f;log 75f;0.30;1f;0.05;0.20);mcCfg));
scenCfg:`forwardShiftPct`spotShiftPct`volShift!(0.05;0.05;0.05);

pos1:`tradeId`commodity`optionSetup`modelInputs`scenarioConfig`quantity`contractMultiplier!(`HW1;`wti;setupCall;inputsAll;scenCfg;1f;1000f);
pos2:@[pos1;(`tradeId;`commodity;`quantity);:;(`HW2;`brent;-1f)];
positions:(pos1;pos2);

greekCfg:.commodity.modelreport.defaultGreekConfig[];
disCfg:.commodity.modelreport.defaultDisagreementConfig[];
limCfg:.commodity.modelreport.defaultModelRiskLimitConfig[];

step1Meta:`runId`runDate`portfolioName!(`HRUN1;2026.05.27;`bookHW);
step1Hist:`limitHistory`summaryHistory!((();()));
step1:.commodity.modelreport.runPortfolioRiskWithLimitHistory[positions;greekCfg;disCfg;limCfg;3;step1Meta;step1Hist];

expectedKeys:`currentRun`limitHistory`summaryHistory`historyDashboard;
.testutil.assertTrue[all expectedKeys in key step1;"step1 has expected keys"];
.testutil.assertTrue[7=count step1`limitHistory;"step1 limitHistory has 7 rows (one snapshot of 7 limits)"];
.testutil.assertTrue[1=count step1`summaryHistory;"step1 summaryHistory has 1 row"];

step2Meta:`runId`runDate`portfolioName!(`HRUN2;2026.05.28;`bookHW);
step2Hist:`limitHistory`summaryHistory!(step1`limitHistory;step1`summaryHistory);
step2:.commodity.modelreport.runPortfolioRiskWithLimitHistory[positions;greekCfg;disCfg;limCfg;3;step2Meta;step2Hist];

.testutil.assertTrue[14=count step2`limitHistory;"step2 limitHistory has 14 rows (2 snapshots)"];
.testutil.assertTrue[2=count step2`summaryHistory;"step2 summaryHistory has 2 rows"];

uniqueRunIds:asc distinct (step2`limitHistory)`runId;
.testutil.assertTrue[uniqueRunIds~asc `HRUN1`HRUN2;"both runIds preserved across history"];

dashboard:step2`historyDashboard;
.testutil.assertTrue[`historySummary in key dashboard;"dashboard has historySummary"];
.testutil.assertTrue[`breachTrend in key dashboard;"dashboard has breachTrend"];
.testutil.assertTrue[dashboard[`dashboardStatus] in `OK`warning`breach`ERROR;"dashboardStatus is a valid status symbol"];

-1 "PASS test_modelreport_limit_history_wrapper: step1Rows=",string[count step1`limitHistory],", step2Rows=",string[count step2`limitHistory];
