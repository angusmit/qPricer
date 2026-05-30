\l core/init.q
/ With loose thresholds, every limit reports OK and overallStatus = OK.
/ Use a real two-position portfolio and override thresholds to be larger than
/ any plausible value.

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

pos1:`tradeId`commodity`optionSetup`modelInputs`scenarioConfig`quantity`contractMultiplier!(`L_OK1;`wti;setupCall;inputsAll;scenCfg;1f;1000f);
pos2:@[pos1;(`tradeId;`commodity;`optionSetup;`quantity);:;(`L_OK2;`brent;setupPut;-2f)];

greekCfg:.commodity.modelreport.defaultGreekConfig[];
disCfg:.commodity.modelreport.defaultDisagreementConfig[];

looseLimits:`grossPriceRangeExposureWarning`grossPriceRangeExposureBreach`maxScenarioPnlRangeWarning`maxScenarioPnlRangeBreach`maxPrimarySensitivityRangeWarning`maxPrimarySensitivityRangeBreach`maxVolatilityVegaRangeWarning`maxVolatilityVegaRangeBreach`maxJumpSensitivityWarning`maxJumpSensitivityBreach`warningAlertCountWarning`warningAlertCountBreach`errorTradeCountWarning`errorTradeCountBreach!(
    1e9;2e9;1e9;2e9;1e9;2e9;1e9;2e9;1e9;2e9;1000;2000;100;200);

portfolioRes:.commodity.modelreport.runPortfolioRisk[(pos1;pos2);greekCfg;disCfg];
limitTbl:.commodity.modelreport.checkPortfolioModelRiskLimits[portfolioRes;looseLimits];
limitSum:.commodity.modelreport.limitStatusSummary limitTbl;

.testutil.assertTrue[7=count limitTbl;"limitTable has seven rows"];
.testutil.assertTrue[all (limitTbl`status)=`OK;"every limit OK under loose thresholds"];
.testutil.assertTrue[`OK=limitSum`overallStatus;"overallStatus OK"];
.testutil.assertTrue[7=limitSum`okCount;"okCount = 7"];
.testutil.assertTrue[0=limitSum`warningCount;"warningCount = 0"];
.testutil.assertTrue[0=limitSum`breachCount;"breachCount = 0"];

breachReport:.commodity.modelreport.limitBreachReport limitTbl;
.testutil.assertTrue[0=count breachReport;"breachReport empty when all OK"];

-1 "PASS test_modelreport_portfolio_limits_ok: okCount=",string[limitSum`okCount];
