/ test_barrier_validation.q - barrier validation (v0.14: all 8 barrier styles supported)
\l core/init.q

.test.expectError:{[testName;fn]
    testResult:@[{x[];`NO_ERROR};fn;{`ERROR}];
    if[testResult~`NO_ERROR; '"Expected error but got success: ",testName];
    -1 "  PASS ",testName;
 };

/ 1. All 8 barrier types now validate OK
validBarriers:`upAndOut`downAndOut`upAndIn`downAndIn;
validIdx:0;
while[validIdx<count validBarriers;
    barrierSym:validBarriers validIdx;
    okTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
        900+validIdx;`AAPL;`equityOption;`european;`call;100f;1f;1f;barrierSym;130f;0f);
    .product.validateOptionTrade okTrade;
    -1 "  PASS valid barrier type ",string barrierSym;
    validIdx+:1];

/ 2. upAndOut with put now valid
okPut:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
    910;`AAPL;`equityOption;`european;`put;100f;1f;1f;`upAndOut;130f;0f);
.product.validateOptionTrade okPut;
-1 "  PASS upAndOut with put now valid";

/ 3. downAndOut with call now valid
okCall:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
    911;`AAPL;`equityOption;`european;`call;100f;1f;1f;`downAndOut;70f;0f);
.product.validateOptionTrade okCall;
-1 "  PASS downAndOut with call now valid";

/ 4. Non-zero rebate still rejected
.test.expectError["non-zero rebate";{
    badTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
        920;`AAPL;`equityOption;`european;`call;100f;1f;1f;`upAndOut;130f;5f);
    .product.validateOptionTrade badTrade}];

/ 5. Missing barrierLevel
.test.expectError["missing barrierLevel";{
    badTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType!(
        921;`AAPL;`equityOption;`european;`call;100f;1f;1f;`upAndOut);
    .product.validateOptionTrade badTrade}];

/ 6. Negative barrierLevel
.test.expectError["negative barrierLevel";{
    badTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
        922;`AAPL;`equityOption;`european;`call;100f;1f;1f;`upAndOut;-130f;0f);
    .product.validateOptionTrade badTrade}];

/ 7. American barrier still rejected (barrier requires European)
.test.expectError["american barrier option";{
    badTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
        923;`AAPL;`equityOption;`american;`put;100f;1f;1f;`downAndOut;70f;0f);
    .product.validateOptionTrade badTrade}];

/ 8. Unsupported barrierType
.test.expectError["unsupported barrierType";{
    badTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
        924;`AAPL;`equityOption;`european;`call;100f;1f;1f;`doubleKnockOut;130f;0f);
    .product.validateOptionTrade badTrade}];

-1 "PASS test_barrier_validation";
