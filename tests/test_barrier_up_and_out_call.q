/ test_barrier_up_and_out_call.q - validate up-and-out call pricing
\l lib/init.q

marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(
    `AAPL;100f;0.05;0f;0.2);
model:.model.createBlackScholesModel[];
config:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;300;4000;0f;300f;`linear;1b;1b);

vanillaCallTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1f);
barrierCallTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
    101;`AAPL;`equityOption;`european;`call;100f;1f;1f;`upAndOut;130f;0f);

vanillaResult:.engine.priceOption[vanillaCallTrade;marketData;model;config];
barrierResult:.engine.priceOption[barrierCallTrade;marketData;model;config];
vanillaUnitPrice:vanillaResult`unitPrice;
barrierUnitPrice:barrierResult`unitPrice;

/ 1. Barrier price >= 0
if[barrierUnitPrice<0f; '"FAIL: barrier call price is negative"];

/ 2. Barrier price < vanilla price
if[not barrierUnitPrice<vanillaUnitPrice; '"FAIL: barrier call should be less than vanilla call"];

/ 3. Materially lower (at least 10% discount)
priceDifference:vanillaUnitPrice-barrierUnitPrice;
if[priceDifference<0.5; '"FAIL: barrier discount too small: ",string priceDifference];

/ 4. Breached barrier: price should be ~0
breachedMkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;131f;0.05;0f;0.2);
breachedResult:.engine.priceOption[barrierCallTrade;breachedMkt;model;config];
breachedUnitPrice:breachedResult`unitPrice;
if[breachedUnitPrice>0.01; '"FAIL: breached barrier should give ~0 price, got ",string breachedUnitPrice];

-1 "PASS test_barrier_up_and_out_call: barrier=",string[barrierUnitPrice],", vanilla=",string[vanillaUnitPrice],", breached=",string breachedUnitPrice;
