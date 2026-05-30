/ test_local_vol_american_put.q - local-vol American put >= local-vol European put
\l core/init.q
lvFn:{[spotValue;timePoint] 0.2};
lvMkt:.market.createLocalVolatilityMarketData[`AAPL;100f;0.05;0f;lvFn];
lvModel:.model.createLocalVolatilityModel[];
cfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);
amPut:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`american;`put;100f;1f;1f);
euPut:@[amPut;`exerciseStyle;:;`european];
amPrice:(.engine.priceOption[amPut;lvMkt;lvModel;cfg])`unitPrice;
euPrice:(.engine.priceOption[euPut;lvMkt;lvModel;cfg])`unitPrice;
if[amPrice<euPrice-0.001; '"FAIL: local-vol American put should be >= European"];
-1 "PASS test_local_vol_american_put: american=",string[amPrice],", european=",string euPrice;
