\l lib/init.q
.test.expectError:{[testName;fn]
    testResult:@[{x[];`NO_ERROR};fn;{`ERROR}];
    if[testResult~`NO_ERROR; '"Expected error but got success: ",testName];
    -1 "  PASS ",testName;
 };
.test.expectError["negative initialVariance";{.bates.validateParams `initialVariance`longRunVariance`meanReversion`volOfVol`correlation`jumpIntensity`jumpMean`jumpVolatility`riskFreeRate`dividendYield!(-0.01;0.04;2.0;0.3;-0.7;0.5;-0.1;0.3;0.05;0f)}];
.test.expectError["zero longRunVariance";{.bates.validateParams `initialVariance`longRunVariance`meanReversion`volOfVol`correlation`jumpIntensity`jumpMean`jumpVolatility`riskFreeRate`dividendYield!(0.04;0.0;2.0;0.3;-0.7;0.5;-0.1;0.3;0.05;0f)}];
.test.expectError["zero meanReversion";{.bates.validateParams `initialVariance`longRunVariance`meanReversion`volOfVol`correlation`jumpIntensity`jumpMean`jumpVolatility`riskFreeRate`dividendYield!(0.04;0.04;0.0;0.3;-0.7;0.5;-0.1;0.3;0.05;0f)}];
.test.expectError["negative volOfVol";{.bates.validateParams `initialVariance`longRunVariance`meanReversion`volOfVol`correlation`jumpIntensity`jumpMean`jumpVolatility`riskFreeRate`dividendYield!(0.04;0.04;2.0;-0.1;-0.7;0.5;-0.1;0.3;0.05;0f)}];
.test.expectError["correlation >= 1";{.bates.validateParams `initialVariance`longRunVariance`meanReversion`volOfVol`correlation`jumpIntensity`jumpMean`jumpVolatility`riskFreeRate`dividendYield!(0.04;0.04;2.0;0.3;1.0;0.5;-0.1;0.3;0.05;0f)}];
.test.expectError["negative jumpIntensity";{.bates.validateParams `initialVariance`longRunVariance`meanReversion`volOfVol`correlation`jumpIntensity`jumpMean`jumpVolatility`riskFreeRate`dividendYield!(0.04;0.04;2.0;0.3;-0.7;-0.5;-0.1;0.3;0.05;0f)}];
.test.expectError["negative jumpVolatility";{.bates.validateParams `initialVariance`longRunVariance`meanReversion`volOfVol`correlation`jumpIntensity`jumpMean`jumpVolatility`riskFreeRate`dividendYield!(0.04;0.04;2.0;0.3;-0.7;0.5;-0.1;-0.3;0.05;0f)}];
-1 "PASS test_bates_validation";
