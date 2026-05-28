/ test_heston_validation.q
\l lib/init.q

.test.expectError:{[testName;fn]
    testResult:@[{x[];`NO_ERROR};fn;{`ERROR}];
    if[testResult~`NO_ERROR; '"Expected error but got success: ",testName];
    -1 "  PASS ",testName;
 };

.test.expectError["negative initialVariance";{
    bad:`initialVariance`longRunVariance`meanReversion`volOfVol`correlation`riskFreeRate`dividendYield!(-0.01;0.04;2.0;0.3;-0.7;0.05;0.0);
    .heston.validateHestonParams bad}];

.test.expectError["zero longRunVariance";{
    bad:`initialVariance`longRunVariance`meanReversion`volOfVol`correlation`riskFreeRate`dividendYield!(0.04;0.0;2.0;0.3;-0.7;0.05;0.0);
    .heston.validateHestonParams bad}];

.test.expectError["zero meanReversion";{
    bad:`initialVariance`longRunVariance`meanReversion`volOfVol`correlation`riskFreeRate`dividendYield!(0.04;0.04;0.0;0.3;-0.7;0.05;0.0);
    .heston.validateHestonParams bad}];

.test.expectError["negative volOfVol";{
    bad:`initialVariance`longRunVariance`meanReversion`volOfVol`correlation`riskFreeRate`dividendYield!(0.04;0.04;2.0;-0.1;-0.7;0.05;0.0);
    .heston.validateHestonParams bad}];

.test.expectError["correlation > 1";{
    bad:`initialVariance`longRunVariance`meanReversion`volOfVol`correlation`riskFreeRate`dividendYield!(0.04;0.04;2.0;0.3;1.5;0.05;0.0);
    .heston.validateHestonParams bad}];

.test.expectError["correlation < -1";{
    bad:`initialVariance`longRunVariance`meanReversion`volOfVol`correlation`riskFreeRate`dividendYield!(0.04;0.04;2.0;0.3;-1.5;0.05;0.0);
    .heston.validateHestonParams bad}];

-1 "PASS test_heston_validation";
