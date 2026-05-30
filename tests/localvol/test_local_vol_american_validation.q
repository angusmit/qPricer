/ test_local_vol_american_validation.q
\l core/init.q

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

/ 1. local-vol American put succeeds
amPut:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`american;`put;100f;1f;1f);
putResult:.engine.priceOption[amPut;lvMkt;lvModel;cfg];
if[not putResult[`unitPrice]>0f; '"FAIL: LV American put should price OK"];
-1 "  PASS local-vol American put prices OK";

/ 2. local-vol American call succeeds
amCall:@[amPut;`optionType;:;`call];
callResult:.engine.priceOption[amCall;lvMkt;lvModel;cfg];
if[not callResult[`unitPrice]>0f; '"FAIL: LV American call should price OK"];
-1 "  PASS local-vol American call prices OK";

/ 3. local-vol + American + barrier still rejected
.test.expectError["local-vol + American + barrier";{
    badTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
        900;`AAPL;`equityOption;`american;`call;100f;1f;1f;`upAndOut;130f;0f);
    .engine.priceOption[badTrade;lvMkt;lvModel;cfg]}];

-1 "PASS test_local_vol_american_validation";
