/ test_local_vol_barrier_down_and_out_put.q
\l core/init.q
lvFn:{[spotValue;timePoint] 0.2+(0f|0.001*(100f-spotValue))};
lvMkt:.market.createLocalVolatilityMarketData[`AAPL;100f;0.05;0f;lvFn];
lvModel:.model.createLocalVolatilityModel[];
cfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;300;4000;0f;300f;`linear;1b;1b);

barrierTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
    1;`AAPL;`equityOption;`european;`put;100f;1f;1f;`downAndOut;70f;0f);
vanillaTrade:@[barrierTrade;`barrierType;:;`none];
barrierPrice:(.engine.priceOption[barrierTrade;lvMkt;lvModel;cfg])`unitPrice;
vanillaPrice:(.engine.priceOption[vanillaTrade;lvMkt;lvModel;cfg])`unitPrice;
if[not barrierPrice>0f; '"FAIL: LV barrier price should be positive"];
if[not barrierPrice<=vanillaPrice+0.001; '"FAIL: LV barrier <= LV vanilla"];
-1 "PASS test_local_vol_barrier_down_and_out_put: barrier=",string[barrierPrice],", vanilla=",string vanillaPrice;
