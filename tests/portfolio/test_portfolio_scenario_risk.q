/ test_portfolio_scenario_risk.q - portfolio scenario risk
\l lib/init.q

tradeTable:([]
    tradeId:1 2 3 4;
    underlying:`AAPL`AAPL`AAPL`AAPL;
    productType:`equityOption`equityOption`equityOption`equityOption;
    exerciseStyle:`european`european`american`european;
    optionType:`call`put`put`call;
    strike:100 100 100 100f;
    expiry:1 1 1 1f;
    notional:1000000 1000000 1000000 1000000f;
    barrierType:`none`none`none`upAndOut;
    barrierLevel:0N 0N 0N 130f;
    rebate:0 0 0 0f);

marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
model:.model.createBlackScholesModel[];
config:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

portfolioScenarioTable:.portfolio.generatePortfolioScenarioReport[tradeTable;marketData;model;config];

/ 1. Row count: 4 trades * 9 scenarios = 36
totalRows:count portfolioScenarioTable;
if[not totalRows=36; '"FAIL: expected 36 rows, got ",string totalRows];

/ 2. All status OK
if[not all portfolioScenarioTable[`status]=`OK; '"FAIL: not all status OK"];

/ 3. Each trade has 9 scenarios
distinctTradeIds:distinct portfolioScenarioTable`tradeId;
if[not 4=count distinctTradeIds; '"FAIL: expected 4 distinct tradeIds"];

/ 4. Base scenario PnL = 0
baseRows:portfolioScenarioTable where portfolioScenarioTable[`scenario]=`base;
if[not 4=count baseRows; '"FAIL: expected 4 base rows"];
if[any (abs baseRows`unitPnL)>0.0001; '"FAIL: base unitPnL should be 0"];
if[any (abs baseRows`notionalPnL)>0.0001; '"FAIL: base notionalPnL should be 0"];

/ 5. Spot up increases European call price
euroCallSpotUp:portfolioScenarioTable where (portfolioScenarioTable[`tradeId]=1) & portfolioScenarioTable[`scenario]=`spotUp1Pct;
euroCallBase:portfolioScenarioTable where (portfolioScenarioTable[`tradeId]=1) & portfolioScenarioTable[`scenario]=`base;
if[not (euroCallSpotUp[`unitPrice]0)>(euroCallBase[`unitPrice]0); '"FAIL: spot up should increase European call"];

/ 6. Spot down decreases European call price
euroCallSpotDn:portfolioScenarioTable where (portfolioScenarioTable[`tradeId]=1) & portfolioScenarioTable[`scenario]=`spotDown1Pct;
if[not (euroCallSpotDn[`unitPrice]0)<(euroCallBase[`unitPrice]0); '"FAIL: spot down should decrease European call"];

scenarioCount:count distinct portfolioScenarioTable`scenario;
-1 "PASS test_portfolio_scenario_risk: rows=",string[totalRows],", scenarios=",string scenarioCount;
