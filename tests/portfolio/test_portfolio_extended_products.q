/ test_portfolio_extended_products.q - portfolio with expanded product set
\l core/init.q

tradeTable:([]
    tradeId:1 2 3 4 5 6;
    underlying:`AAPL`AAPL`AAPL`AAPL`AAPL`AAPL;
    productType:6#`equityOption;
    exerciseStyle:`european`european`american`american`european`european;
    optionType:`call`put`put`call`call`put;
    strike:100 100 100 100 100 100f;
    expiry:1 1 1 1 1 1f;
    notional:1000000 1000000 1000000 1000000 1000000 1000000f;
    barrierType:`none`none`none`none`upAndOut`downAndIn;
    barrierLevel:0N 0N 0N 0N 130 70f;
    rebate:0 0 0 0 0 0f);

mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0.01;0.2);
mdl:.model.createBlackScholesModel[];
cfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

pricingResult:.portfolio.priceTrades[tradeTable;mkt;mdl;cfg];
okCount:sum pricingResult[`status]=`OK;
if[not okCount=6; '"FAIL: all 6 trades should price OK, got ",string okCount];
if[any (pricingResult`unitPrice)<=0f; '"FAIL: all prices should be positive"];
-1 "PASS test_portfolio_extended_products: rows=6, okRows=",string okCount;
