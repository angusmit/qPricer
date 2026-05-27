/ test_sabr_validation.q
\l lib/init.q

.test.expectError:{[testName;fn]
    testResult:@[{x[];`NO_ERROR};fn;{`ERROR}];
    if[testResult~`NO_ERROR; '"Expected error but got success: ",testName];
    -1 "  PASS ",testName;
 };

.test.expectError["negative alpha";{.sabr.validateParams `alpha`beta`rho`nu!(-0.1;0.5;-0.3;0.4)}];
.test.expectError["beta < 0";{.sabr.validateParams `alpha`beta`rho`nu!(0.2;-0.1;-0.3;0.4)}];
.test.expectError["beta > 1";{.sabr.validateParams `alpha`beta`rho`nu!(0.2;1.5;-0.3;0.4)}];
.test.expectError["rho <= -1";{.sabr.validateParams `alpha`beta`rho`nu!(0.2;0.5;-1.0;0.4)}];
.test.expectError["rho >= 1";{.sabr.validateParams `alpha`beta`rho`nu!(0.2;0.5;1.0;0.4)}];
.test.expectError["negative nu";{.sabr.validateParams `alpha`beta`rho`nu!(0.2;0.5;-0.3;-0.1)}];
.test.expectError["zero forward";{.sabr.impliedVol[0f;100f;1f;`alpha`beta`rho`nu!(0.2;0.5;-0.3;0.4)]}];
.test.expectError["zero strike";{.sabr.impliedVol[100f;0f;1f;`alpha`beta`rho`nu!(0.2;0.5;-0.3;0.4)]}];
.test.expectError["zero expiry";{.sabr.impliedVol[100f;100f;0f;`alpha`beta`rho`nu!(0.2;0.5;-0.3;0.4)]}];
.test.expectError["empty smile";{.sabr.validateMarketSmile ()}];
.test.expectError["missing marketIV";{.sabr.validateMarketSmile enlist `strike`expiry!(100f;1f)}];

-1 "PASS test_sabr_validation";
