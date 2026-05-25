/ test_local_vol_validation.q - local vol scope and error handling
\l lib/init.q

.test.expectError:{[testName;fn]
    testResult:@[{x[];`NO_ERROR};fn;{`ERROR}];
    if[testResult~`NO_ERROR; '"Expected error but got success: ",testName];
    -1 "  PASS ",testName;
 };

baseTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1f);
localVolFunction:{[spotValue;timePoint] 0.2};
localVolMkt:.market.createLocalVolatilityMarketData[`AAPL;100f;0.05;0f;localVolFunction];
lvModel:.model.createLocalVolatilityModel[];
explicitCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);
cnCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;200;500;0f;300f;`linear;1b;1b);

/ 1. Local vol with Crank-Nicolson should error
.test.expectError["local vol with CN";{.engine.priceOption[baseTrade;localVolMkt;lvModel;cnCfg]}];

/ 2. Local vol with American put should error
americanPut:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    900;`AAPL;`equityOption;`american;`put;100f;1f;1f);
.test.expectError["local vol with American";{.engine.priceOption[americanPut;localVolMkt;lvModel;explicitCfg]}];

/ 3. Local vol with barrier should error
barrierTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
    901;`AAPL;`equityOption;`european;`call;100f;1f;1f;`upAndOut;130f;0f);
.test.expectError["local vol with barrier";{.engine.priceOption[barrierTrade;localVolMkt;lvModel;explicitCfg]}];

/ 4. Missing local vol function should error
.test.expectError["missing local vol function";{
    badMkt:`underlying`spot`riskFreeRate`dividendYield!(`AAPL;100f;0.05;0f);
    .market.validateLocalVolatilityMarketData badMkt}];

/ 5. Local vol returning zero should error
.test.expectError["local vol returning zero";{
    zeroVolFn:{[spotValue;timePoint] 0f};
    .market.createLocalVolatilityMarketData[`AAPL;100f;0.05;0f;zeroVolFn]}];

/ 6. Local vol returning negative should error
.test.expectError["local vol returning negative";{
    negVolFn:{[spotValue;timePoint] -0.1};
    .market.createLocalVolatilityMarketData[`AAPL;100f;0.05;0f;negVolFn]}];

-1 "PASS test_local_vol_validation";
