\l lib/init.q
/ runEnsemble returns one summary row per (strategy, path).

pathCfg:`spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.25;5;1f%252f;0.02;0f;42);
ensemble:.strategy.path.ensemble[pathCfg;3;100];
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `E1;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;80;150;0f;300f;`linear;1b;1b);

gsCfg:.strategy.defaultConfig `gammaScalp;
gsCfg:@[gsCfg;`stepYears;:;1f%252f];
svCfg:.strategy.defaultConfig `shortVariance;
svCfg:@[svCfg;(`forecastVol;`entryMargin;`stepYears);:;(0.10;0.02;1f%252f)];
specs:((`gammaScalp;gsCfg);(`shortVariance;svCfg));

summaryTbl:.strategy.portfolio.runEnsemble[specs;trade;ensemble;bsModel;fdmCfg];

requiredCols:`strategyName`pathId`totalPnl`maxDrawdown`numRebalances`status`errorMessage;
.testutil.assertTableColumns[summaryTbl;requiredCols;"ensemble summary schema"];
.testutil.assertTrue[6=count summaryTbl;"2 strategies x 3 paths = 6 rows"];
.testutil.assertTrue[3=count summaryTbl where (summaryTbl`strategyName)=`gammaScalp;"3 gammaScalp rows"];
.testutil.assertTrue[3=count summaryTbl where (summaryTbl`strategyName)=`shortVariance;"3 shortVariance rows"];
.testutil.assertTrue[(asc 0 1 2)~asc distinct summaryTbl`pathId;"pathIds 0..2"];
.testutil.assertTrue[all (summaryTbl`status)=`OK;"every cell OK"];

emptyEnsResult:.[.strategy.portfolio.runEnsemble;(specs;trade;();bsModel;fdmCfg);{`ERROR}];
.testutil.assertTrue[emptyEnsResult~`ERROR;"empty ensemble rejected"];

-1 "PASS test_strategy_portfolio_run_ensemble: rows=",string[count summaryTbl];
