\l lib/init.q
/ portfolioDisagreementExposure aggregates per-position disagreement metrics
/ into gross/net/max measures. Verify counts and non-negativity of gross/max
/ metrics on a small portfolio.

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

pos1:`tradeId`commodity`optionSetup`modelInputs`scenarioConfig`quantity`contractMultiplier!(`P_A;`wti;setupCall;inputsAll;scenCfg;1f;1000f);
pos2:@[pos1;(`tradeId;`commodity;`optionSetup;`quantity);:;(`P_B;`brent;setupPut;-2f)];

positions:(pos1;pos2);
greekCfg:.commodity.modelreport.defaultGreekConfig[];
disCfg:.commodity.modelreport.defaultDisagreementConfig[];

portfolioResult:.commodity.modelreport.runPortfolioRisk[positions;greekCfg;disCfg];
exposure:.commodity.modelreport.portfolioDisagreementExposure portfolioResult;

expectedKeys:`tradeCount`okTradeCount`warningTradeCount`errorTradeCount`grossPriceRangeExposure`netScenarioPnlRange`grossScenarioPnlRange`maxPriceRange`maxPrimarySensitivityRange`maxVolatilityVegaRange`maxScenarioPnlRange`maxJumpSensitivity`status`errorMessage;
.testutil.assertTrue[all expectedKeys in key exposure;"exposure has expected keys"];

.testutil.assertTrue[2=exposure`tradeCount;"tradeCount = 2"];
.testutil.assertTrue[2=exposure`okTradeCount;"okTradeCount = 2"];
.testutil.assertTrue[0=exposure`errorTradeCount;"errorTradeCount = 0"];
.testutil.assertTrue[exposure[`warningTradeCount]>=0;"warningTradeCount non-negative"];
.testutil.assertTrue[exposure[`warningTradeCount]<=2;"warningTradeCount <= tradeCount"];
.testutil.assertTrue[`OK=exposure`status;"status OK"];

.testutil.assertTrue[exposure[`grossPriceRangeExposure]>=0f;"grossPriceRangeExposure non-negative"];
.testutil.assertTrue[exposure[`grossScenarioPnlRange]>=0f;"grossScenarioPnlRange non-negative"];
.testutil.assertTrue[exposure[`maxPriceRange]>=0f;"maxPriceRange non-negative"];
.testutil.assertTrue[exposure[`maxPrimarySensitivityRange]>=0f;"maxPrimarySensitivityRange non-negative"];
.testutil.assertTrue[exposure[`maxVolatilityVegaRange]>=0f;"maxVolatilityVegaRange non-negative"];
.testutil.assertTrue[exposure[`maxScenarioPnlRange]>=0f;"maxScenarioPnlRange non-negative"];

posSummaryTbl:portfolioResult`positionSummary;
expectedGross:sum (abs posSummaryTbl`priceRange)*(abs posSummaryTbl`quantity)*abs posSummaryTbl`contractMultiplier;
.testutil.assertNear[exposure`grossPriceRangeExposure;expectedGross;1e-9;"grossPriceRangeExposure = sum |range|*|qty|*|mult|"];

expectedMaxPrice:max abs posSummaryTbl`priceRange;
.testutil.assertNear[exposure`maxPriceRange;expectedMaxPrice;1e-12;"maxPriceRange = max |priceRange|"];

-1 "PASS test_modelreport_portfolio_exposure: gross=",string[exposure`grossPriceRangeExposure],", maxScenPnl=",string[exposure`maxScenarioPnlRange];
