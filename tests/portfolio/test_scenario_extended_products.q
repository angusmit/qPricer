/ test_scenario_extended_products.q - scenario risk for expanded product set
\l lib/init.q
tradeTable:([]
    tradeId:1 2 3 4;
    underlying:`AAPL`AAPL`AAPL`AAPL;
    productType:4#`equityOption;
    exerciseStyle:`european`american`american`european;
    optionType:`call`put`call`call;
    strike:100 100 100 100f;
    expiry:1 1 1 1f;
    notional:1000000 1000000 1000000 1000000f;
    barrierType:`none`none`none`upAndOut;
    barrierLevel:0N 0N 0N 130f;
    rebate:0 0 0 0f);

mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0.01;0.2);
mdl:.model.createBlackScholesModel[];
cfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

scenarioResult:.portfolio.generatePortfolioScenarioReport[tradeTable;mkt;mdl;cfg];
expectedRows:4*9;
actualRows:count scenarioResult;
if[not actualRows=expectedRows; '"FAIL: expected ",string[expectedRows]," rows, got ",string actualRows];
-1 "PASS test_scenario_extended_products: rows=",string actualRows;
