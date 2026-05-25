/ test_barrier_validation.q - barrier-specific input validation
\l lib/init.q

.test.expectError:{[testName;fn]
    testResult:@[{x[];`NO_ERROR};fn;{`ERROR}];
    if[testResult~`NO_ERROR; '"Expected error but got success: ",testName];
    -1 "  PASS ",testName;
 };

baseMkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
baseModel:.model.createBlackScholesModel[];
baseCfg:.config.createFiniteDifferenceConfig[`method`numberOfSpotSteps`numberOfTimeSteps`maximumSpot!(`explicit;200;2000;300f)];

/ 1. Unsupported barrierType
.test.expectError["unsupported barrierType upAndIn";{
    badTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
        900;`AAPL;`equityOption;`european;`call;100f;1f;1f;`upAndIn;130f;0f);
    .product.validateOptionTrade badTrade}];

/ 2. upAndOut with put
.test.expectError["upAndOut with put";{
    badTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
        901;`AAPL;`equityOption;`european;`put;100f;1f;1f;`upAndOut;130f;0f);
    .product.validateOptionTrade badTrade}];

/ 3. downAndOut with call
.test.expectError["downAndOut with call";{
    badTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
        902;`AAPL;`equityOption;`european;`call;100f;1f;1f;`downAndOut;70f;0f);
    .product.validateOptionTrade badTrade}];

/ 4. Non-zero rebate
.test.expectError["non-zero rebate";{
    badTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
        903;`AAPL;`equityOption;`european;`call;100f;1f;1f;`upAndOut;130f;5f);
    .product.validateOptionTrade badTrade}];

/ 5. Missing barrierLevel
.test.expectError["missing barrierLevel";{
    badTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType!(
        904;`AAPL;`equityOption;`european;`call;100f;1f;1f;`upAndOut);
    .product.validateOptionTrade badTrade}];

/ 6. Negative barrierLevel
.test.expectError["negative barrierLevel";{
    badTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
        905;`AAPL;`equityOption;`european;`call;100f;1f;1f;`upAndOut;-130f;0f);
    .product.validateOptionTrade badTrade}];

/ 7. American barrier (not supported in v0.5)
.test.expectError["american barrier option";{
    badTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
        906;`AAPL;`equityOption;`american;`put;100f;1f;1f;`downAndOut;70f;0f);
    .product.validateOptionTrade badTrade}];

-1 "PASS test_barrier_validation";
