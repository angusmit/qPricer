/ test_local_vol_american_call.q - local-vol American call >= local-vol European call
\l lib/init.q
lvFn:{[spotValue;timePoint] 0.2};
lvMkt:.market.createLocalVolatilityMarketData[`AAPL;100f;0.05;0.03;lvFn];
lvModel:.model.createLocalVolatilityModel[];
cfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);
amCall:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`american;`call;100f;1f;1f);
euCall:@[amCall;`exerciseStyle;:;`european];
amPrice:(.engine.priceOption[amCall;lvMkt;lvModel;cfg])`unitPrice;
euPrice:(.engine.priceOption[euCall;lvMkt;lvModel;cfg])`unitPrice;
if[amPrice<euPrice-0.001; '"FAIL: local-vol American call should be >= European"];
-1 "PASS test_local_vol_american_call: american=",string[amPrice],", european=",string euPrice;
