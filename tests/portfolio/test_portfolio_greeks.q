/ test_portfolio_greeks.q - portfolio Greeks (v0.14: all products get Greeks)
\l core/init.q

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

portfolioGreeksTable:.portfolio.calculatePortfolioGreeks[tradeTable;marketData;model;config];

if[not 4=count portfolioGreeksTable; '"FAIL: expected 4 rows"];

/ v0.14: all 4 trades should produce OK Greeks
statusValues:portfolioGreeksTable`status;
okRowCount:sum statusValues=`OK;
if[not okRowCount=4; '"FAIL: expected all 4 rows OK, got ",string okRowCount];

/ European call delta > 0
if[not (portfolioGreeksTable`delta)0 > 0f; '"FAIL: European call delta > 0"];
/ European put delta < 0
if[not (portfolioGreeksTable`delta)1 < 0f; '"FAIL: European put delta < 0"];
/ American put delta < 0
if[not (portfolioGreeksTable`delta)2 < 0f; '"FAIL: American put delta < 0"];
/ Barrier call delta > 0
if[not (portfolioGreeksTable`delta)3 > 0f; '"FAIL: barrier call delta > 0"];

-1 "PASS test_portfolio_greeks: rows=4, okRows=",string okRowCount;
