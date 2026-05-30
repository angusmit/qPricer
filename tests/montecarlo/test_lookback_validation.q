/ test_lookback_validation.q
\l core/init.q

.test.expectError:{[testName;fn]
    testResult:@[{x[];`NO_ERROR};fn;{`ERROR}];
    if[testResult~`NO_ERROR; '"Expected error but got success: ",testName];
    -1 "  PASS ",testName;
 };

.test.expectError["invalid lookbackStyle";{
    bad:`tradeId`underlying`productType`exerciseStyle`optionType`lookbackStyle`strike`expiry`notional`observationCount!(
        1;`AAPL;`lookbackOption;`european;`call;`barrier;100f;1f;1f;50);
    .lookback.validateLookbackTrade bad}];

.test.expectError["American lookback";{
    bad:`tradeId`underlying`productType`exerciseStyle`optionType`lookbackStyle`strike`expiry`notional`observationCount!(
        1;`AAPL;`lookbackOption;`american;`call;`fixed;100f;1f;1f;50);
    .lookback.validateLookbackTrade bad}];

.test.expectError["negative strike fixed";{
    bad:`tradeId`underlying`productType`exerciseStyle`optionType`lookbackStyle`strike`expiry`notional`observationCount!(
        1;`AAPL;`lookbackOption;`european;`call;`fixed;-100f;1f;1f;50);
    .lookback.validateLookbackTrade bad}];

.test.expectError["zero expiry";{
    bad:`tradeId`underlying`productType`exerciseStyle`optionType`lookbackStyle`strike`expiry`notional`observationCount!(
        1;`AAPL;`lookbackOption;`european;`call;`fixed;100f;0f;1f;50);
    .lookback.validateLookbackTrade bad}];

.test.expectError["invalid observationCount";{
    bad:`tradeId`underlying`productType`exerciseStyle`optionType`lookbackStyle`strike`expiry`notional`observationCount!(
        1;`AAPL;`lookbackOption;`european;`call;`fixed;100f;1f;1f;0);
    .lookback.validateLookbackTrade bad}];

-1 "PASS test_lookback_validation";
