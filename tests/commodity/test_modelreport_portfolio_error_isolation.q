\l lib/init.q
/ A bad position (invalid optionType inside its optionSetup) must not crash
/ the portfolio. The good position should still produce all its rows; the bad
/ position should appear with positionStatus=`ERROR and an explanatory
/ errorMessage; portfolioSummary.status should be `warning.

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

goodPos:`tradeId`commodity`optionSetup`modelInputs`scenarioConfig`quantity`contractMultiplier!(`OK_POS;`wti;setupCall;inputsAll;scenCfg;1f;1000f);
badPos:@[goodPos;(`tradeId;`commodity;`optionSetup);:;(`BAD_POS;`gas;setupBad)];

positions:(goodPos;badPos);
greekCfg:.commodity.modelreport.defaultGreekConfig[];
disCfg:.commodity.modelreport.defaultDisagreementConfig[];

portfolioResult:.commodity.modelreport.runPortfolioRisk[positions;greekCfg;disCfg];

posSummaryTbl:portfolioResult`positionSummary;
.testutil.assertTrue[2=count posSummaryTbl;"positionSummary has both positions"];

okRowsSelect:posSummaryTbl where (posSummaryTbl`positionStatus)=`OK;
.testutil.assertTrue[1=count okRowsSelect;"one OK position"];
.testutil.assertTrue[`OK_POS=first okRowsSelect`tradeId;"OK position tradeId preserved"];

errRowsSelect:posSummaryTbl where (posSummaryTbl`positionStatus)=`ERROR;
.testutil.assertTrue[1=count errRowsSelect;"one ERROR position"];
.testutil.assertTrue[`BAD_POS=first errRowsSelect`tradeId;"bad position tradeId preserved"];
.testutil.assertTrue[0<count first errRowsSelect`positionErrorMessage;"bad position has non-empty errorMessage"];

portfolioSummary:portfolioResult`portfolioSummary;
.testutil.assertTrue[2=portfolioSummary`tradeCount;"tradeCount = 2"];
.testutil.assertTrue[1=portfolioSummary`okTradeCount;"okTradeCount = 1"];
.testutil.assertTrue[1=portfolioSummary`errorTradeCount;"errorTradeCount = 1"];
.testutil.assertTrue[`warning=portfolioSummary`status;"portfolioSummary status warning (mixed)"];

priceRowsTbl:portfolioResult`priceRows;
.testutil.assertTrue[4=count priceRowsTbl;"only OK position contributes 4 price rows"];
.testutil.assertTrue[all (priceRowsTbl`tradeId)=`OK_POS;"all priceRows tagged OK_POS"];

alertRowsTbl:portfolioResult`alertRows;
.testutil.assertTrue[6=count alertRowsTbl;"only OK position contributes 6 alert rows"];

exposure:.commodity.modelreport.portfolioDisagreementExposure portfolioResult;
.testutil.assertTrue[1=exposure`okTradeCount;"exposure okTradeCount = 1"];
.testutil.assertTrue[1=exposure`errorTradeCount;"exposure errorTradeCount = 1"];

-1 "PASS test_modelreport_portfolio_error_isolation: status=",string[portfolioSummary`status],", okTradeCount=",string[portfolioSummary`okTradeCount];
