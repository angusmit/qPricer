\l lib/init.q
/ runPortfolioRiskWithLimits stitches runPortfolioDashboard with the limit
/ checker. Verify dictionary keys, expected limit names in limitTable, and
/ that limitSummary exposes overallStatus.

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

pos1:`tradeId`commodity`optionSetup`modelInputs`scenarioConfig`quantity`contractMultiplier!(`RW1;`wti;setupCall;inputsAll;scenCfg;1f;1000f);
pos2:@[pos1;(`tradeId;`commodity;`optionSetup;`quantity);:;(`RW2;`brent;setupPut;-2f)];

greekCfg:.commodity.modelreport.defaultGreekConfig[];
disCfg:.commodity.modelreport.defaultDisagreementConfig[];
limCfg:.commodity.modelreport.defaultModelRiskLimitConfig[];

fullResult:.commodity.modelreport.runPortfolioRiskWithLimits[(pos1;pos2);greekCfg;disCfg;limCfg;2];

expectedKeys:`portfolioRisk`dashboard`limitTable`limitSummary`breachReport;
.testutil.assertTrue[all expectedKeys in key fullResult;"runPortfolioRiskWithLimits has all expected keys"];

limitTbl:fullResult`limitTable;
expectedLimits:`grossPriceRangeExposure`maxScenarioPnlRange`maxPrimarySensitivityRange`maxVolatilityVegaRange`maxJumpSensitivity`warningAlertCount`errorTradeCount;
.testutil.assertTrue[all expectedLimits in limitTbl`limitName;"limitTable contains all expected limit names"];
.testutil.assertTrue[7=count limitTbl;"limitTable has seven rows"];

limitSum:fullResult`limitSummary;
.testutil.assertTrue[7=limitSum`limitCount;"limitSummary.limitCount = 7"];
.testutil.assertTrue[limitSum[`overallStatus] in `OK`warning`breach`ERROR;"overallStatus is a valid status symbol"];

breachReport:fullResult`breachReport;
.testutil.assertTrue[98h=type breachReport;"breachReport is a table"];

dashboardKeys:`portfolioSummary`alertSummary`disagreementExposure`worstOffenders;
.testutil.assertTrue[all dashboardKeys in key fullResult`dashboard;"dashboard preserved inside wrapper"];

portfolioRisk:fullResult`portfolioRisk;
.testutil.assertTrue[`positionSummary in key portfolioRisk;"portfolioRisk includes positionSummary"];

-1 "PASS test_modelreport_portfolio_risk_with_limits: limitCount=",string[limitSum`limitCount],", overall=",string[limitSum`overallStatus];
