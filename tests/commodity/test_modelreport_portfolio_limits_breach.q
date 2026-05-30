\l core/init.q
/ With tight breach thresholds, at least one limit should breach.
/ Verify overallStatus = breach, breachAmount = metric - breachThreshold,
/ and that limitBreachReport surfaces the breach rows.

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

pos1:`tradeId`commodity`optionSetup`modelInputs`scenarioConfig`quantity`contractMultiplier!(`L_B1;`wti;setupCall;inputsAll;scenCfg;1f;1000f);
pos2:@[pos1;(`tradeId;`commodity;`optionSetup;`quantity);:;(`L_B2;`brent;setupPut;-2f)];

greekCfg:.commodity.modelreport.defaultGreekConfig[];
disCfg:.commodity.modelreport.defaultDisagreementConfig[];

tightBreachLimits:`grossPriceRangeExposureWarning`grossPriceRangeExposureBreach`maxScenarioPnlRangeWarning`maxScenarioPnlRangeBreach`maxPrimarySensitivityRangeWarning`maxPrimarySensitivityRangeBreach`maxVolatilityVegaRangeWarning`maxVolatilityVegaRangeBreach`maxJumpSensitivityWarning`maxJumpSensitivityBreach`warningAlertCountWarning`warningAlertCountBreach`errorTradeCountWarning`errorTradeCountBreach!(
    0.001;0.01;0.001;0.01;0.001;0.01;0.001;0.01;0.001;0.01;0;1;0;1);

portfolioRes:.commodity.modelreport.runPortfolioRisk[(pos1;pos2);greekCfg;disCfg];
limitTbl:.commodity.modelreport.checkPortfolioModelRiskLimits[portfolioRes;tightBreachLimits];
limitSum:.commodity.modelreport.limitStatusSummary limitTbl;

.testutil.assertTrue[limitSum[`breachCount]>0;"at least one breach under tight breach thresholds"];
.testutil.assertTrue[`breach=limitSum`overallStatus;"overallStatus = breach"];

breachRows:limitTbl where (limitTbl`status)=`breach;
.testutil.assertTrue[all (breachRows`breachAmount)>0f;"breach rows have positive breachAmount"];

priceExpRow:limitTbl first where (limitTbl`limitName)=`grossPriceRangeExposure;
.testutil.assertNear[priceExpRow`breachAmount;priceExpRow[`metricValue]-priceExpRow`breachThreshold;1e-9;"breachAmount = metric - breachThreshold"];

breachReport:.commodity.modelreport.limitBreachReport limitTbl;
.testutil.assertTrue[(count breachReport)>=limitSum`breachCount;"breachReport contains all breach rows"];
firstReportRow:breachReport 0;
.testutil.assertTrue[`breach=firstReportRow`status;"breachReport first row is breach (sorted ERROR/breach/warning)"];

-1 "PASS test_modelreport_portfolio_limits_breach: breachCount=",string[limitSum`breachCount],", overall=",string[limitSum`overallStatus];
