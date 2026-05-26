/ test_local_vol_barrier_validation.q
\l lib/init.q

.test.expectError:{[testName;fn]
    testResult:@[{x[];`NO_ERROR};fn;{`ERROR}];
    if[testResult~`NO_ERROR; '"Expected error but got success: ",testName];
    -1 "  PASS ",testName;
 };

lvFn:{[spotValue;timePoint] 0.2};
lvMkt:.market.createLocalVolatilityMarketData[`AAPL;100f;0.05;0f;lvFn];
lvModel:.model.createLocalVolatilityModel[];
cfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

/ 1. LV + European barrier succeeds
euBarrier:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1f;`upAndOut;130f;0f);
okResult:.engine.priceOption[euBarrier;lvMkt;lvModel;cfg];
if[not okResult[`unitPrice]>0f; '"FAIL: LV European barrier should price OK"];
-1 "  PASS LV + European barrier prices OK";

/ 2. LV + American + barrier rejected (barrier requires European at product level)
.test.expectError["LV + American + barrier";{
    badTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
        900;`AAPL;`equityOption;`american;`call;100f;1f;1f;`upAndOut;130f;0f);
    .engine.priceOption[badTrade;lvMkt;lvModel;cfg]}];

-1 "PASS test_local_vol_barrier_validation";
