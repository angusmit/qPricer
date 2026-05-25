/ test_implied_vol_invalid_prices.q - IV solver rejects invalid inputs
\l lib/init.q

.test.expectError:{[testName;fn]
    functionResult:@[{x[];`NO_ERROR};fn;{`ERROR}];
    if[functionResult~`NO_ERROR; '"Expected error but got success: ",testName];
    -1 "  PASS ",testName;
 };

/ 1. Negative market price
.test.expectError["negative market price";
    {.iv.impliedVolatility[`call;-1f;100f;100f;1f;0.05;0f]}];

/ 2. Zero market price
.test.expectError["zero market price";
    {.iv.impliedVolatility[`call;0f;100f;100f;1f;0.05;0f]}];

/ 3. Price below no-arbitrage lower bound
.test.expectError["price below lower bound";
    {.iv.impliedVolatility[`call;0.0001;100f;100f;1f;0.05;0f]}];

/ 4. Price above no-arbitrage upper bound
.test.expectError["price above upper bound";
    {.iv.impliedVolatility[`call;200f;100f;100f;1f;0.05;0f]}];

/ 5. Unsupported option type
.test.expectError["unsupported option type";
    {.iv.impliedVolatility[`digital;5f;100f;100f;1f;0.05;0f]}];

/ 6. Bad solver config: non-positive lower bound
.test.expectError["non-positive lower vol bound";
    {badCfg:`lowerVolatilityBound`upperVolatilityBound`tolerance`maximumIterations!(0f;5f;1e-8;100);
     .iv.impliedVolatilityWithConfig[`call;10f;100f;100f;1f;0.05;0f;badCfg]}];

/ 7. Bad solver config: upper <= lower
.test.expectError["upper vol bound <= lower";
    {badCfg:`lowerVolatilityBound`upperVolatilityBound`tolerance`maximumIterations!(0.5;0.1;1e-8;100);
     .iv.impliedVolatilityWithConfig[`call;10f;100f;100f;1f;0.05;0f;badCfg]}];

/ 8. Bad solver config: non-positive tolerance
.test.expectError["non-positive tolerance";
    {badCfg:`lowerVolatilityBound`upperVolatilityBound`tolerance`maximumIterations!(0.01;5f;0f;100);
     .iv.impliedVolatilityWithConfig[`call;10f;100f;100f;1f;0.05;0f;badCfg]}];

-1 "PASS test_implied_vol_invalid_prices";
