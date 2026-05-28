\l lib/init.q
/ validatePortfolioPositions accepts well-formed portfolios and rejects
/ positions that are missing required keys or that pass non-numeric quantity
/ / contractMultiplier.

mcCfg:.montecarlo.defaultMcConfig[];
mcCfg:@[mcCfg;`pathCount;:;1000];
mcCfg:@[mcCfg;`timeStepCount;:;10];

setupCall:`optionType`strikePrice`expiry`riskFreeRate`forwardPrice`spotPrice!(`call;75f;1f;0.05;75f;75f);
inputsAll:`black76`schwartz`schwartz2`mrjump!(
    enlist[`volatility]!enlist 0.30;
    `x0`params!(log 75f;`meanReversionSpeed`longRunLogMean`volatility!(1.5;log 75f;0.30));
    `shortFactor0`longFactor0`params!(0f;log 75f;`meanReversionSpeed`shortVolatility`longVolatility`longDrift`correlation!(1.5;0.30;0.15;0f;0.25));
    `x0`params`mcConfig!(log 75f;`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility!(2f;log 75f;0.30;1f;0.05;0.20);mcCfg));
scenCfg:`forwardShiftPct`spotShiftPct`volShift!(0.05;0.05;0.05);

okPos1:`tradeId`commodity`optionSetup`modelInputs`scenarioConfig`quantity`contractMultiplier!(`T1;`wti;setupCall;inputsAll;scenCfg;1f;1000f);
okPos2:@[okPos1;`tradeId;:;`T2];
.commodity.modelreport.validatePortfolioPositions[(okPos1;okPos2)];
-1 "  validate passes on two valid positions";

emptyResult:@[.commodity.modelreport.validatePortfolioPositions;();{`ERROR}];
.testutil.assertTrue[emptyResult~`ERROR;"empty portfolio rejected"];

missingTradeId:(okPos1;(`commodity`optionSetup`modelInputs`scenarioConfig`quantity`contractMultiplier)#okPos2);
missingResult:@[.commodity.modelreport.validatePortfolioPositions;missingTradeId;{`ERROR}];
.testutil.assertTrue[missingResult~`ERROR;"missing tradeId rejected"];

badQty:@[okPos1;`quantity;:;`badQty];
badQtyResult:@[.commodity.modelreport.validatePortfolioPositions;enlist badQty;{`ERROR}];
.testutil.assertTrue[badQtyResult~`ERROR;"non-numeric quantity rejected"];

badMult:@[okPos2;`contractMultiplier;:;"oops"];
badMultResult:@[.commodity.modelreport.validatePortfolioPositions;enlist badMult;{`ERROR}];
.testutil.assertTrue[badMultResult~`ERROR;"non-numeric contractMultiplier rejected"];

-1 "PASS test_modelreport_portfolio_positions";
