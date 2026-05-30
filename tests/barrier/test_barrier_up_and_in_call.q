/ test_barrier_up_and_in_call.q
\l core/init.q
knockInTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1f;`upAndIn;130f;0f);
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
mdl:.model.createBlackScholesModel[];
cfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;300;4000;0f;300f;`linear;1b;1b);
knockInPrice:(.engine.priceOption[knockInTrade;mkt;mdl;cfg])`unitPrice;
if[not knockInPrice>=0f; '"FAIL: knock-in price should be >= 0"];
-1 "PASS test_barrier_up_and_in_call: price=",string knockInPrice;
