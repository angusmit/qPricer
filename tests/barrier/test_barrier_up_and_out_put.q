/ test_barrier_up_and_out_put.q
\l lib/init.q
barrierTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
    1;`AAPL;`equityOption;`european;`put;100f;1f;1f;`upAndOut;130f;0f);
vanillaTrade:@[barrierTrade;`barrierType;:;`none];
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
mdl:.model.createBlackScholesModel[];
cfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;300;4000;0f;300f;`linear;1b;1b);
barrierPrice:(.engine.priceOption[barrierTrade;mkt;mdl;cfg])`unitPrice;
vanillaPrice:(.engine.priceOption[vanillaTrade;mkt;mdl;cfg])`unitPrice;
if[not barrierPrice>0f; '"FAIL: barrier price should be positive"];
if[not barrierPrice<=vanillaPrice+0.001; '"FAIL: barrier <= vanilla"];
-1 "PASS test_barrier_up_and_out_put: barrier=",string[barrierPrice],", vanilla=",string vanillaPrice;
