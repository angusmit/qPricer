\l core/init.q
/ With tight warning thresholds but very loose breach thresholds, every metric
/ should report warning and none should breach. Verify overallStatus warning
/ and that breachAmount is 0 for warning rows.

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

pos1:`tradeId`commodity`optionSetup`modelInputs`scenarioConfig`quantity`contractMultiplier!(`L_W1;`wti;setupCall;inputsAll;scenCfg;1f;1000f);
pos2:@[pos1;(`tradeId;`commodity;`optionSetup;`quantity);:;(`L_W2;`brent;setupPut;-2f)];

greekCfg:.commodity.modelreport.defaultGreekConfig[];
disCfg:.commodity.modelreport.defaultDisagreementConfig[];

tightWarningLimits:`grossPriceRangeExposureWarning`grossPriceRangeExposureBreach`maxScenarioPnlRangeWarning`maxScenarioPnlRangeBreach`maxPrimarySensitivityRangeWarning`maxPrimarySensitivityRangeBreach`maxVolatilityVegaRangeWarning`maxVolatilityVegaRangeBreach`maxJumpSensitivityWarning`maxJumpSensitivityBreach`warningAlertCountWarning`warningAlertCountBreach`errorTradeCountWarning`errorTradeCountBreach!(
    0.01;1e9;0.01;1e9;0.01;1e9;0.01;1e9;0.01;1e9;0;1000;0;200);

portfolioRes:.commodity.modelreport.runPortfolioRisk[(pos1;pos2);greekCfg;disCfg];
limitTbl:.commodity.modelreport.checkPortfolioModelRiskLimits[portfolioRes;tightWarningLimits];
limitSum:.commodity.modelreport.limitStatusSummary limitTbl;

.testutil.assertTrue[0=limitSum`breachCount;"zero breaches under loose breach thresholds"];
.testutil.assertTrue[limitSum[`warningCount]>0;"at least one warning under tight warning thresholds"];
.testutil.assertTrue[`warning=limitSum`overallStatus;"overallStatus = warning"];

warningRows:limitTbl where (limitTbl`status)=`warning;
.testutil.assertTrue[all 0f=warningRows`breachAmount;"warning rows have breachAmount = 0"];

breachReport:.commodity.modelreport.limitBreachReport limitTbl;
.testutil.assertTrue[(count breachReport)=limitSum`warningCount;"breachReport row count matches warningCount when no breaches"];
.testutil.assertTrue[all `warning=breachReport`status;"breachReport rows all warning"];

-1 "PASS test_modelreport_portfolio_limits_warning: warningCount=",string[limitSum`warningCount];
