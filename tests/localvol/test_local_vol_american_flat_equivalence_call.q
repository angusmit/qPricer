/ test_local_vol_american_flat_equivalence_call.q
\l core/init.q
lvFn:{[spotValue;timePoint] 0.2};
lvMkt:.market.createLocalVolatilityMarketData[`AAPL;100f;0.05;0.03;lvFn];
flatMkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0.03;0.2);
lvModel:.model.createLocalVolatilityModel[];
bsModel:.model.createBlackScholesModel[];
cfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);
amCall:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`american;`call;100f;1f;1f);
lvPrice:(.engine.priceOption[amCall;lvMkt;lvModel;cfg])`unitPrice;
flatPrice:(.engine.priceOption[amCall;flatMkt;bsModel;cfg])`unitPrice;
priceDiff:abs lvPrice-flatPrice;
if[priceDiff>0.001; '"FAIL: flat local-vol should match constant vol, diff=",string priceDiff];
-1 "PASS test_local_vol_american_flat_equivalence_call: diff=",string priceDiff;
