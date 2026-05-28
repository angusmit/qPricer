/ test_local_vol_barrier_flat_equivalence.q - flat LV barrier = constant vol barrier
\l lib/init.q
lvFn:{[spotValue;timePoint] 0.2};
lvMkt:.market.createLocalVolatilityMarketData[`AAPL;100f;0.05;0f;lvFn];
flatMkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
lvModel:.model.createLocalVolatilityModel[];
bsModel:.model.createBlackScholesModel[];
cfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;300;4000;0f;300f;`linear;1b;1b);

uoCall:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1f;`upAndOut;130f;0f);
lvPrice:(.engine.priceOption[uoCall;lvMkt;lvModel;cfg])`unitPrice;
bsPrice:(.engine.priceOption[uoCall;flatMkt;bsModel;cfg])`unitPrice;
priceDiff:abs lvPrice-bsPrice;
if[priceDiff>0.001; '"FAIL: flat LV barrier should match BS barrier, diff=",string priceDiff];
-1 "PASS test_local_vol_barrier_flat_equivalence: diff=",string priceDiff;
