/ test_local_vol_validation.q - local vol scope and error handling (v0.14)
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

/ 1. Local vol with CN still errors
.test.expectError["local vol with CN";{.engine.priceOption[baseTrade;localVolMkt;lvModel;cnCfg]}];

/ 2. Local vol with American now succeeds (v0.14)
americanPut:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    900;`AAPL;`equityOption;`american;`put;100f;1f;1f);
lvAmericanResult:.engine.priceOption[americanPut;localVolMkt;lvModel;explicitCfg];
if[not lvAmericanResult[`unitPrice]>0f; '"FAIL: local vol American put should price"];
-1 "  PASS local vol with American now supported";

/ 3. Local vol with barrier now succeeds (v0.14)
barrierTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
    901;`AAPL;`equityOption;`european;`call;100f;1f;1f;`upAndOut;130f;0f);
lvBarrierResult:.engine.priceOption[barrierTrade;localVolMkt;lvModel;explicitCfg];
if[not lvBarrierResult[`unitPrice]>0f; '"FAIL: local vol barrier should price"];
-1 "  PASS local vol with barrier now supported";

/ 4. Local vol + American + barrier still errors
.test.expectError["local vol + American + barrier";{
    lvAmBarrier:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
        902;`AAPL;`equityOption;`american;`call;100f;1f;1f;`upAndOut;130f;0f);
    .engine.priceOption[lvAmBarrier;localVolMkt;lvModel;explicitCfg]}];

/ 5. Missing local vol function still errors
.test.expectError["missing local vol function";{
    badMkt:`underlying`spot`riskFreeRate`dividendYield!(`AAPL;100f;0.05;0f);
    .market.validateLocalVolatilityMarketData badMkt}];

/ 6. Local vol returning zero still errors
.test.expectError["local vol returning zero";{
    zeroVolFn:{[spotValue;timePoint] 0f};
    .market.createLocalVolatilityMarketData[`AAPL;100f;0.05;0f;zeroVolFn]}];

/ 7. Local vol returning negative still errors
.test.expectError["local vol returning negative";{
    negVolFn:{[spotValue;timePoint] -0.1};
    .market.createLocalVolatilityMarketData[`AAPL;100f;0.05;0f;negVolFn]}];

-1 "PASS test_local_vol_validation";
