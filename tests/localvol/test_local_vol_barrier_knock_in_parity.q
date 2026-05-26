/ test_local_vol_barrier_knock_in_parity.q - KI + KO = vanilla under flat LV
\l lib/init.q
lvFn:{[spotValue;timePoint] 0.2};
lvMkt:.market.createLocalVolatilityMarketData[`AAPL;100f;0.05;0f;lvFn];
lvModel:.model.createLocalVolatilityModel[];
cfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;300;4000;0f;300f;`linear;1b;1b);

vanillaCall:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1f;`none;0Nf;0f);
uoCall:@[@[vanillaCall;`barrierType;:;`upAndOut];`barrierLevel;:;130f];
uiCall:@[@[vanillaCall;`barrierType;:;`upAndIn];`barrierLevel;:;130f];

vPrice:(.engine.priceOption[vanillaCall;lvMkt;lvModel;cfg])`unitPrice;
koPrice:(.engine.priceOption[uoCall;lvMkt;lvModel;cfg])`unitPrice;
kiPrice:(.engine.priceOption[uiCall;lvMkt;lvModel;cfg])`unitPrice;
parityError:abs (kiPrice+koPrice)-vPrice;
if[parityError>0.01; '"FAIL: LV barrier parity error=",string parityError];
-1 "PASS test_local_vol_barrier_knock_in_parity: vanilla=",string[vPrice],", KO+KI=",string kiPrice+koPrice;
