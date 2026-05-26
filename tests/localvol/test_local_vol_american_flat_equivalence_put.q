/ test_local_vol_american_flat_equivalence_put.q
\l lib/init.q
lvFn:{[spotValue;timePoint] 0.2};
lvMkt:.market.createLocalVolatilityMarketData[`AAPL;100f;0.05;0f;lvFn];
flatMkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
lvModel:.model.createLocalVolatilityModel[];
bsModel:.model.createBlackScholesModel[];
cfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);
amPut:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`american;`put;100f;1f;1f);
lvPrice:(.engine.priceOption[amPut;lvMkt;lvModel;cfg])`unitPrice;
flatPrice:(.engine.priceOption[amPut;flatMkt;bsModel;cfg])`unitPrice;
priceDiff:abs lvPrice-flatPrice;
if[priceDiff>0.001; '"FAIL: flat local-vol should match constant vol, diff=",string priceDiff];
-1 "PASS test_local_vol_american_flat_equivalence_put: diff=",string priceDiff;
