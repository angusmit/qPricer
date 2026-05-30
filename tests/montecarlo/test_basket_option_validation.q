/ test_basket_option_validation.q
\l core/init.q

.test.expectError:{[testName;fn]
    testResult:@[{x[];`NO_ERROR};fn;{`ERROR}];
    if[testResult~`NO_ERROR; '"Expected error but got success: ",testName];
    -1 "  PASS ",testName;
 };

/ 1. Mismatched symbols/weights
.test.expectError["mismatched symbols/weights";{
    bad:`tradeId`productType`basketSymbols`basketWeights`optionType`exerciseStyle`strike`expiry`notional!(
        1;`basketOption;`AAPL`MSFT;enlist 0.5;`call;`european;100f;1f;1f);
    .basket.validateBasketTrade bad}];

/ 2. Empty basket
.test.expectError["empty basket";{
    bad:`tradeId`productType`basketSymbols`basketWeights`optionType`exerciseStyle`strike`expiry`notional!(
        1;`basketOption;`symbol$();`float$();`call;`european;100f;1f;1f);
    .basket.validateBasketTrade bad}];

/ 3. Invalid optionType
.test.expectError["invalid optionType";{
    bad:`tradeId`productType`basketSymbols`basketWeights`optionType`exerciseStyle`strike`expiry`notional!(
        1;`basketOption;`AAPL`MSFT;0.5 0.5;`digital;`european;100f;1f;1f);
    .basket.validateBasketTrade bad}];

/ 4. American basket rejected
.test.expectError["American basket";{
    bad:`tradeId`productType`basketSymbols`basketWeights`optionType`exerciseStyle`strike`expiry`notional!(
        1;`basketOption;`AAPL`MSFT;0.5 0.5;`call;`american;100f;1f;1f);
    .basket.validateBasketTrade bad}];

/ 5. Negative strike
.test.expectError["negative strike";{
    bad:`tradeId`productType`basketSymbols`basketWeights`optionType`exerciseStyle`strike`expiry`notional!(
        1;`basketOption;`AAPL`MSFT;0.5 0.5;`call;`european;-100f;1f;1f);
    .basket.validateBasketTrade bad}];

-1 "PASS test_basket_option_validation";
