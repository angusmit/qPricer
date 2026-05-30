\l core/init.q
/ A portfolio with one bad position still produces a limit table. With tight
/ errorTradeCount thresholds the failed position should drive the
/ errorTradeCount limit to warning or breach without crashing.

mcCfg:.montecarlo.defaultMcConfig[];
mcCfg:@[mcCfg;`pathCount;:;3000];
mcCfg:@[mcCfg;`timeStepCount;:;15];

setupCall:`optionType`strikePrice`expiry`riskFreeRate`forwardPrice`spotPrice!(`call;75f;1f;0.05;75f;75f);
setupBad:@[setupCall;`optionType;:;`straddle];
inputsAll:`black76`schwartz`schwartz2`mrjump!(
    enlist[`volatility]!enlist 0.30;
    `x0`params!(log 75f;`meanReversionSpeed`longRunLogMean`volatility!(1.5;log 75f;0.30));
    `shortFactor0`longFactor0`params!(0f;log 75f;`meanReversionSpeed`shortVolatility`longVolatility`longDrift`correlation!(1.5;0.30;0.15;0f;0.25));
    `x0`params`mcConfig!(log 75f;`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility!(2f;log 75f;0.30;1f;0.05;0.20);mcCfg));
scenCfg:`forwardShiftPct`spotShiftPct`volShift!(0.05;0.05;0.05);

goodPos:`tradeId`commodity`optionSetup`modelInputs`scenarioConfig`quantity`contractMultiplier!(`L_GOOD;`wti;setupCall;inputsAll;scenCfg;1f;1000f);
badPos:@[goodPos;(`tradeId;`commodity;`optionSetup);:;(`L_BAD;`gas;setupBad)];

greekCfg:.commodity.modelreport.defaultGreekConfig[];
disCfg:.commodity.modelreport.defaultDisagreementConfig[];

tightErrorLimits:`grossPriceRangeExposureWarning`grossPriceRangeExposureBreach`maxScenarioPnlRangeWarning`maxScenarioPnlRangeBreach`maxPrimarySensitivityRangeWarning`maxPrimarySensitivityRangeBreach`maxVolatilityVegaRangeWarning`maxVolatilityVegaRangeBreach`maxJumpSensitivityWarning`maxJumpSensitivityBreach`warningAlertCountWarning`warningAlertCountBreach`errorTradeCountWarning`errorTradeCountBreach!(
    1e9;2e9;1e9;2e9;1e9;2e9;1e9;2e9;1e9;2e9;1000;2000;1;2);

portfolioRes:.commodity.modelreport.runPortfolioRisk[(goodPos;badPos);greekCfg;disCfg];

.testutil.assertTrue[`warning=portfolioRes[`portfolioSummary;`status];"portfolio status warning with mixed positions"];
.testutil.assertTrue[1=portfolioRes[`portfolioSummary;`errorTradeCount];"one ERROR trade"];

limitTbl:.commodity.modelreport.checkPortfolioModelRiskLimits[portfolioRes;tightErrorLimits];
.testutil.assertTrue[7=count limitTbl;"limitTable has seven rows even with failed position"];

errorTradeRow:limitTbl first where (limitTbl`limitName)=`errorTradeCount;
.testutil.assertTrue[1f=errorTradeRow`metricValue;"errorTradeCount metricValue = 1"];
.testutil.assertTrue[errorTradeRow[`status] in `warning`breach;"errorTradeCount limit triggers warning or breach"];

limitSum:.commodity.modelreport.limitStatusSummary limitTbl;
.testutil.assertTrue[limitSum[`overallStatus] in `warning`breach;"overallStatus warning or breach when error count tripped"];

-1 "PASS test_modelreport_portfolio_limits_error_isolation: errorRow status=",string[errorTradeRow`status],", overall=",string[limitSum`overallStatus];
