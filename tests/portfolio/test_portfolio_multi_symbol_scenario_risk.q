/ test_portfolio_multi_symbol_scenario_risk.q - multi-symbol scenario risk
\l lib/init.q

tradeTable:([]
    tradeId:1 2 3;
    underlying:`AAPL`MSFT`NVDA;
    productType:`equityOption`equityOption`equityOption;
    exerciseStyle:`european`european`european;
    optionType:`call`put`call;
    strike:100 250 800f;
    expiry:1 1 1f;
    notional:1000000 1000000 1000000f;
    barrierType:`none`none`none;
    barrierLevel:0N 0N 0N;
    rebate:0 0 0f);

spotTable:([] underlying:`AAPL`MSFT`NVDA; spot:100 250 800f);
volatilityTable:([] underlying:`AAPL`MSFT`NVDA; volatility:0.20 0.25 0.35);
rateTable:([] expiry:enlist 1f; riskFreeRate:enlist 0.05);
dividendTable:([] underlying:`AAPL`MSFT`NVDA; dividendYield:0 0.01 0f);

marketDataBook:.marketbook.createMarketDataBook[spotTable;volatilityTable;rateTable;dividendTable];
model:.model.createBlackScholesModel[];
config:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;200;500;0f;2400f;`linear;1b;1b);

portfolioScenarioTable:.portfolio.generatePortfolioScenarioReportWithMarketDataBook[tradeTable;marketDataBook;model;config];

totalRows:count portfolioScenarioTable;
if[not totalRows=27; '"FAIL: expected 27 rows, got ",string totalRows];
if[not all portfolioScenarioTable[`status]=`OK; '"FAIL: not all status OK"];

baseRows:portfolioScenarioTable where portfolioScenarioTable[`scenario]=`base;
if[not 3=count baseRows; '"FAIL: expected 3 base rows"];
if[any (abs baseRows`unitPnL)>0.0001; '"FAIL: base unitPnL should be 0"];

aaplSpotUp:portfolioScenarioTable where (portfolioScenarioTable[`tradeId]=1) & portfolioScenarioTable[`scenario]=`spotUp1Pct;
aaplBase:portfolioScenarioTable where (portfolioScenarioTable[`tradeId]=1) & portfolioScenarioTable[`scenario]=`base;
if[not (aaplSpotUp[`unitPrice]0)>(aaplBase[`unitPrice]0); '"FAIL: AAPL call spot up should increase price"];

msftSpotUp:portfolioScenarioTable where (portfolioScenarioTable[`tradeId]=2) & portfolioScenarioTable[`scenario]=`spotUp1Pct;
msftBase:portfolioScenarioTable where (portfolioScenarioTable[`tradeId]=2) & portfolioScenarioTable[`scenario]=`base;
if[not (msftSpotUp[`unitPrice]0)<(msftBase[`unitPrice]0); '"FAIL: MSFT put spot up should decrease price"];

scenarioCount:count distinct portfolioScenarioTable`scenario;
-1 "PASS test_portfolio_multi_symbol_scenario_risk: rows=",string[totalRows],", scenarios=",string scenarioCount;
