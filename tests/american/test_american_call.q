/ test_american_call.q - American call basic pricing
\l lib/init.q
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`american;`call;100f;1f;1f);
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0.03;0.2);
mdl:.model.createBlackScholesModel[];
cfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);
priceResult:.engine.priceOption[trade;mkt;mdl;cfg];
amCallPrice:priceResult`unitPrice;
/ American call with dividend should be positive
if[not amCallPrice>0f; '"FAIL: American call price should be positive"];
/ American call should be >= intrinsic
intrinsicVal:0f|mkt[`spot]-trade`strike;
if[amCallPrice<intrinsicVal-0.001; '"FAIL: American call price should be >= intrinsic"];
-1 "PASS test_american_call: price=",string amCallPrice;
