\l core/init.q
/ portfolioRiskDashboard composes portfolioSummary + alertSummary +
/ disagreementExposure + worstOffenders. runPortfolioDashboard is the
/ end-to-end wrapper that also runs the portfolio.

mcCfg:.montecarlo.defaultMcConfig[];
mcCfg:@[mcCfg;`pathCount;:;3000];
mcCfg:@[mcCfg;`timeStepCount;:;15];

setupCall:`optionType`strikePrice`expiry`riskFreeRate`forwardPrice`spotPrice!(`call;75f;1f;0.05;75f;75f);
setupPut:@[setupCall;(`optionType;`strikePrice);:;(`put;80f)];
inputsAll:`black76`schwartz`schwartz2`mrjump!(
    enlist[`volatility]!enlist 0.30;
    `x0`params!(log 75f;`meanReversionSpeed`longRunLogMean`volatility!(1.5;log 75f;0.30));
    `shortFactor0`longFactor0`params!(0f;log 75f;`meanReversionSpeed`shortVolatility`longVolatility`longDrift`correlation!(1.5;0.30;0.15;0f;0.25));
    `x0`params`mcConfig!(log 75f;`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility!(2f;log 75f;0.30;1f;0.05;0.20);mcCfg));
scenCfg:`forwardShiftPct`spotShiftPct`volShift!(0.05;0.05;0.05);

pos1:`tradeId`commodity`optionSetup`modelInputs`scenarioConfig`quantity`contractMultiplier!(`D1;`wti;setupCall;inputsAll;scenCfg;1f;1000f);
pos2:@[pos1;(`tradeId;`commodity;`optionSetup;`quantity);:;(`D2;`brent;setupPut;-2f)];

positions:(pos1;pos2);
greekCfg:.commodity.modelreport.defaultGreekConfig[];
disCfg:.commodity.modelreport.defaultDisagreementConfig[];

portfolioResult:.commodity.modelreport.runPortfolioRisk[positions;greekCfg;disCfg];
dashboard:.commodity.modelreport.portfolioRiskDashboard portfolioResult;

expectedKeys:`portfolioSummary`alertSummary`disagreementExposure`worstOffenders;
.testutil.assertTrue[all expectedKeys in key dashboard;"portfolioRiskDashboard has expected keys"];

portfolioSummary:dashboard`portfolioSummary;
.testutil.assertTrue[2=portfolioSummary`tradeCount;"dashboard portfolioSummary tradeCount = 2"];

alertSummary:dashboard`alertSummary;
.testutil.assertTrue[98h=type alertSummary;"alertSummary is a table"];
.testutil.assertTrue[0<count alertSummary;"alertSummary has rows when alerts present"];

exposure:dashboard`disagreementExposure;
.testutil.assertTrue[`OK=exposure`status;"exposure status OK"];

offenders:dashboard`worstOffenders;
.testutil.assertTrue[0<count offenders;"offenders table populated with topN=3 default"];

runDashboard:.commodity.modelreport.runPortfolioDashboard[positions;greekCfg;disCfg;2];
runDashboardKeys:`portfolioRisk`portfolioSummary`alertSummary`disagreementExposure`worstOffenders;
.testutil.assertTrue[all runDashboardKeys in key runDashboard;"runPortfolioDashboard has expected keys"];

runOffenders:runDashboard`worstOffenders;
metricCount:count distinct runOffenders`metricName;
.testutil.assertTrue[(2*metricCount)=count runOffenders;"runPortfolioDashboard topN=2 -> 2 rows per metric"];

-1 "PASS test_modelreport_portfolio_dashboard: dashboardKeys=",string[count key dashboard],", offenderRows=",string[count offenders];
