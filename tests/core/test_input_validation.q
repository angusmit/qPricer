/ test_input_validation.q - edge-case input validation tests
\l lib/init.q

.test.expectError:{[testName;fn]
    testResult:@[{x[];`NO_ERROR};fn;{`ERROR}];
    if[testResult~`NO_ERROR; '"Expected error but got success: ",testName];
    -1 "  PASS ",testName;
 };

/ --- Base inputs (all valid) ---

baseTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1f);
baseMkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(
    `AAPL;100f;0.05;0f;0.2);
baseModel:.model.createBlackScholesModel[];
baseCfg:.config.createFiniteDifferenceConfig[`method`numberOfSpotSteps`numberOfTimeSteps`maximumSpot!(
    `explicit;200;2000;300f)];

/ --- Trade validation ---
.test.expectError["negative strike";{.product.validateOptionTrade @[baseTrade;`strike;:;-100f]}];
.test.expectError["zero strike";{.product.validateOptionTrade @[baseTrade;`strike;:;0f]}];
.test.expectError["zero expiry";{.product.validateOptionTrade @[baseTrade;`expiry;:;0f]}];
.test.expectError["negative expiry";{.product.validateOptionTrade @[baseTrade;`expiry;:;-1f]}];
.test.expectError["unsupported optionType";{.product.validateOptionTrade @[baseTrade;`optionType;:;`digital]}];
.test.expectError["unsupported exerciseStyle";{.product.validateOptionTrade @[baseTrade;`exerciseStyle;:;`bermudan]}];
.test.expectError["unsupported productType";{.product.validateOptionTrade @[baseTrade;`productType;:;`fxOption]}];
.test.expectError["zero notional";{.product.validateOptionTrade @[baseTrade;`notional;:;0f]}];
.test.expectError["missing strike field";{.product.validateOptionTrade `strike _ baseTrade}];

/ --- Market data validation ---
.test.expectError["negative spot";{.market.validateFlatMarketData @[baseMkt;`spot;:;-100f]}];
.test.expectError["zero volatility";{.market.validateFlatMarketData @[baseMkt;`volatility;:;0f]}];
.test.expectError["negative volatility";{.market.validateFlatMarketData @[baseMkt;`volatility;:;-0.2]}];

/ --- Config validation ---
.test.expectError["unsupported method";{.config.validateFiniteDifferenceConfig @[baseCfg;`method;:;`implicit]}];
.test.expectError["unsupported interpolationMethod";{.config.validateFiniteDifferenceConfig @[baseCfg;`interpolationMethod;:;`cubic]}];
.test.expectError["zero spotSteps";{.config.validateFiniteDifferenceConfig @[baseCfg;`numberOfSpotSteps;:;0]}];
.test.expectError["zero timeSteps";{.config.validateFiniteDifferenceConfig @[baseCfg;`numberOfTimeSteps;:;0]}];
.test.expectError["maximumSpot <= minimumSpot";{.config.validateFiniteDifferenceConfig `method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(`explicit;200;2000;300f;100f;`linear;1b;1b)}];

/ --- Engine validation ---
.test.expectError["maximumSpot below spot";{.engine.priceOption[baseTrade;baseMkt;baseModel;@[baseCfg;`maximumSpot;:;50f]]}];
.test.expectError["maximumSpot below strike";{.engine.priceOption[@[baseTrade;`strike;:;400f];baseMkt;baseModel;baseCfg]}];
.test.expectError["underlying mismatch";{.engine.priceOption[baseTrade;@[baseMkt;`underlying;:;`MSFT];baseModel;baseCfg]}];
.test.expectError["spot outside grid range";{.engine.priceOption[baseTrade;@[baseMkt;`spot;:;400f];baseModel;baseCfg]}];

/ --- CN-specific validation ---
cnCfg:.config.createFiniteDifferenceConfig[`method`numberOfSpotSteps`numberOfTimeSteps`maximumSpot!(`crankNicolson;200;500;300f)];
americanTrade:@[baseTrade;`exerciseStyle;:;`american];
americanTrade:@[americanTrade;`optionType;:;`put];
.test.expectError["CN with American put";{.engine.priceOption[americanTrade;baseMkt;baseModel;cnCfg]}];

barrierTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
    900;`AAPL;`equityOption;`european;`call;100f;1f;1f;`upAndOut;130f;0f);
.test.expectError["CN with barrier option";{.engine.priceOption[barrierTrade;baseMkt;baseModel;cnCfg]}];

-1 "PASS test_input_validation";
