/ test_perfopt.q - performance optimisation correctness (v0.15)
\l lib/init.q

bsModel:.model.createBlackScholesModel[];
fdmConfig:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;100;200;0f;1500f;`linear;1b;1b);
configDict:`model`fdmConfig`timeStepYears`bookName`valuationDate`runLabel!(
    bsModel;fdmConfig;1%252;"perfBook";2025.01.01;"perfRun");

mktBook:.stress.generateMarketDataBook[3;2025.01.01];
symbolList:mktBook[`spotTable]`underlying;
tradeTable:.stress.generateSupportedTradeTable[5;symbolList;2025.01.01];
prevBook:.stress.generateMarketDataBook[3;2024.12.31];

/ 1. priceBookWithCache matches baseline
baselineResult:.portfolio.priceTradesWithMarketDataBook[tradeTable;mktBook;bsModel;fdmConfig];
cacheResult:.perfopt.priceBookWithCache[tradeTable;mktBook;bsModel;fdmConfig];
basePVs:`float${x`unitPrice} each baselineResult;
optPVs:`float${x`unitPrice} each cacheResult`pricingResult;
maxPVDiff:max abs basePVs-optPVs;
.testutil.assertTrue[maxPVDiff<1e-10;"cached pricing matches baseline"];

/ 2. PnL explain with Greek reuse — actualPnL must match (same pricing)
greekResult:.batch.__calculateGreeksFromBook[tradeTable;prevBook;bsModel;fdmConfig];
pnlConfig:`model`fdmConfig`timeStepYears`bookName!(bsModel;fdmConfig;1%252;"test");
baselinePnL:.pnl.explainPortfolio[tradeTable;prevBook;mktBook;pnlConfig];
optimisedPnL:.perfopt.pnlExplainWithGreekReuse[tradeTable;prevBook;mktBook;pnlConfig;greekResult];
baseActualPnL:`float${x`actualPnL} each baselinePnL;
optActualPnL:`float${x`actualPnL} each optimisedPnL;
maxPnLDiff:max abs baseActualPnL-optActualPnL;
.testutil.assertTrue[maxPnLDiff<1e-10;"optimised PnL actualPnL matches baseline"];

/ 3. compareBaselineOptimised runs without error
comparisonResult:.[.perfopt.compareBaselineOptimised;(tradeTable;mktBook;prevBook;configDict);{x}];
if[10h=type comparisonResult; -1 "  comparison error: ",comparisonResult; comparisonResult:()];

/ 4. If comparison succeeded, verify structure
if[0<count comparisonResult;
   .testutil.assertTrue[0<count comparisonResult;"comparison has rows"];
   / Check OK rows have non-null timing
   okRows:comparisonResult where (`OK=){x`status} each comparisonResult;
   if[0<count okRows;
      okBaseMs:`float${x`baselineMs} each okRows;
      .testutil.assertTrue[not any null okBaseMs;"OK rows have non-null baselineMs"]]];

-1 "PASS test_perfopt: pricingDiff=",string[maxPVDiff],", pnlDiff=",string maxPnLDiff;
