\l core/init.q
/ runPortfolioRisk runs runPositionRisk across positions and concatenates the
/ tagged tables. Row counts should scale with the number of positions and the
/ portfolioSummary should report OK / OK trade counts when all positions price.

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
pos3:@[pos1;(`tradeId;`commodity;`quantity);:;(`P_C;`gas;3f)];

positions:(pos1;pos2;pos3);
greekCfg:.commodity.modelreport.defaultGreekConfig[];
disCfg:.commodity.modelreport.defaultDisagreementConfig[];

portfolioResult:.commodity.modelreport.runPortfolioRisk[positions;greekCfg;disCfg];

expectedKeys:`positionResults`priceRows`greeksRows`scenarioPnlRows`alertRows`positionSummary`portfolioSummary`status`errorMessage;
.testutil.assertTrue[all expectedKeys in key portfolioResult;"runPortfolioRisk has expected keys"];

.testutil.assertTrue[3=count portfolioResult`positionResults;"three position results"];

priceRowsTbl:portfolioResult`priceRows;
greeksRowsTbl:portfolioResult`greeksRows;
scenarioPnlRowsTbl:portfolioResult`scenarioPnlRows;
alertRowsTbl:portfolioResult`alertRows;
positionSummaryTbl:portfolioResult`positionSummary;

.testutil.assertTrue[12=count priceRowsTbl;"priceRows = 3 positions x 4 models"];
.testutil.assertTrue[12=count greeksRowsTbl;"greeksRows = 3 positions x 4 models"];
.testutil.assertTrue[12=count scenarioPnlRowsTbl;"scenarioPnlRows = 3 positions x 4 models"];
.testutil.assertTrue[18=count alertRowsTbl;"alertRows = 3 positions x 6 alert types"];
.testutil.assertTrue[3=count positionSummaryTbl;"positionSummary one row per position"];

uniqueTradesInPrice:asc distinct priceRowsTbl`tradeId;
.testutil.assertTrue[uniqueTradesInPrice~asc `P_A`P_B`P_C;"priceRows cover all tradeIds"];
uniqueTradesInAlerts:asc distinct alertRowsTbl`tradeId;
.testutil.assertTrue[uniqueTradesInAlerts~asc `P_A`P_B`P_C;"alertRows cover all tradeIds"];

portfolioSummary:portfolioResult`portfolioSummary;
.testutil.assertTrue[3=portfolioSummary`tradeCount;"portfolioSummary tradeCount = 3"];
.testutil.assertTrue[3=portfolioSummary`okTradeCount;"portfolioSummary okTradeCount = 3"];
.testutil.assertTrue[0=portfolioSummary`errorTradeCount;"portfolioSummary errorTradeCount = 0"];
.testutil.assertTrue[`OK=portfolioSummary`status;"portfolioSummary status OK"];
.testutil.assertTrue[`OK=portfolioResult`status;"portfolio top-level status OK"];

-1 "PASS test_modelreport_run_portfolio_risk: positions=",string[count positions],", priceRows=",string[count priceRowsTbl],", alertRows=",string[count alertRowsTbl];
