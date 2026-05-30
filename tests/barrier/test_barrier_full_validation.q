/ test_barrier_full_validation.q - extended barrier validation
\l core/init.q

.test.expectError:{[testName;fn]
    testResult:@[{x[];`NO_ERROR};fn;{`ERROR}];
    if[testResult~`NO_ERROR; '"Expected error but got success: ",testName];
    -1 "  PASS ",testName;
 };

mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
mdl:.model.createBlackScholesModel[];
cfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

/ 1. All 8 barrier/option combos price without error
combos:((`upAndOut;`call;130f);(`upAndOut;`put;130f);(`downAndOut;`call;70f);(`downAndOut;`put;70f);
        (`upAndIn;`call;130f);(`upAndIn;`put;130f);(`downAndIn;`call;70f);(`downAndIn;`put;70f));
comboIdx:0;
while[comboIdx<count combos;
    combo:combos comboIdx;
    testTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
        100+comboIdx;`AAPL;`equityOption;`european;combo 1;100f;1f;1f;combo 0;combo 2;0f);
    testResult:.[.engine.priceOption;(testTrade;mkt;mdl;cfg);{x}];
    if[10h=type testResult; '"FAIL: ",string[combo 0]," ",string[combo 1]," failed: ",testResult];
    -1 "  PASS ",string[combo 0]," ",string combo 1;
    comboIdx+:1];

/ 2. Invalid barrier type still fails
.test.expectError["invalid barrier type";{
    bad:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
        200;`AAPL;`equityOption;`european;`call;100f;1f;1f;`tripleBarrier;130f;0f);
    .product.validateOptionTrade bad}];

-1 "PASS test_barrier_full_validation";
