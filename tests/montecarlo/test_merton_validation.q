/ test_merton_validation.q
\l lib/init.q
.test.expectError:{[testName;fn]
    testResult:@[{x[];`NO_ERROR};fn;{`ERROR}];
    if[testResult~`NO_ERROR; '"Expected error but got success: ",testName];
    -1 "  PASS ",testName;
 };

.test.expectError["negative volatility";{.merton.validateParams `volatility`jumpIntensity`jumpMean`jumpVolatility`riskFreeRate`dividendYield!(-0.2;0.5;-0.1;0.3;0.05;0f)}];
.test.expectError["negative jumpIntensity";{.merton.validateParams `volatility`jumpIntensity`jumpMean`jumpVolatility`riskFreeRate`dividendYield!(0.2;-0.5;-0.1;0.3;0.05;0f)}];
.test.expectError["negative jumpVolatility";{.merton.validateParams `volatility`jumpIntensity`jumpMean`jumpVolatility`riskFreeRate`dividendYield!(0.2;0.5;-0.1;-0.3;0.05;0f)}];

-1 "PASS test_merton_validation";
