/ test_barrier_down_and_out_put.q - validate down-and-out put pricing
\l lib/init.q

marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(
    `AAPL;100f;0.05;0f;0.2);
model:.model.createBlackScholesModel[];
config:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;300;4000;0f;300f;`linear;1b;1b);

vanillaPutTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    2;`AAPL;`equityOption;`european;`put;100f;1f;1f);
barrierPutTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
    102;`AAPL;`equityOption;`european;`put;100f;1f;1f;`downAndOut;70f;0f);

vanillaResult:.engine.priceOption[vanillaPutTrade;marketData;model;config];
barrierResult:.engine.priceOption[barrierPutTrade;marketData;model;config];
vanillaUnitPrice:vanillaResult`unitPrice;
barrierUnitPrice:barrierResult`unitPrice;

/ 1. Barrier price >= 0
if[barrierUnitPrice<0f; '"FAIL: barrier put price is negative"];

/ 2. Barrier price < vanilla price
if[not barrierUnitPrice<vanillaUnitPrice; '"FAIL: barrier put should be less than vanilla put"];

/ 3. Materially lower
priceDifference:vanillaUnitPrice-barrierUnitPrice;
if[priceDifference<0.1; '"FAIL: barrier discount too small: ",string priceDifference];

/ 4. Breached barrier: price should be ~0
breachedMkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;69f;0.05;0f;0.2);
breachedResult:.engine.priceOption[barrierPutTrade;breachedMkt;model;config];
breachedUnitPrice:breachedResult`unitPrice;
if[breachedUnitPrice>0.01; '"FAIL: breached barrier should give ~0 price, got ",string breachedUnitPrice];

-1 "PASS test_barrier_down_and_out_put: barrier=",string[barrierUnitPrice],", vanilla=",string[vanillaUnitPrice],", breached=",string breachedUnitPrice;
