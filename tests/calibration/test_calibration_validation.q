/ test_calibration_validation.q
\l lib/init.q

.test.expectError:{[testName;fn]
    testResult:@[{x[];`NO_ERROR};fn;{`ERROR}];
    if[testResult~`NO_ERROR; '"Expected error but got success: ",testName];
    -1 "  PASS ",testName;
 };

/ Empty table
.test.expectError["empty table";{
    .calibration.validateMarketOptionTable ()}];

/ Missing column
.test.expectError["missing column";{
    bad:enlist `optionType`strike`expiry!(`;100f;1f);
    .calibration.validateMarketOptionTable bad}];

/ Negative strike
.test.expectError["negative strike";{
    bad:enlist `optionType`strike`expiry`marketPrice!(`call;-100f;1f;10f);
    .calibration.validateMarketOptionTable bad}];

/ Zero expiry
.test.expectError["zero expiry";{
    bad:enlist `optionType`strike`expiry`marketPrice!(`call;100f;0f;10f);
    .calibration.validateMarketOptionTable bad}];

/ Negative marketPrice
.test.expectError["negative marketPrice";{
    bad:enlist `optionType`strike`expiry`marketPrice!(`call;100f;1f;-10f);
    .calibration.validateMarketOptionTable bad}];

/ Invalid optionType
.test.expectError["invalid optionType";{
    bad:enlist `optionType`strike`expiry`marketPrice!(`digital;100f;1f;10f);
    .calibration.validateMarketOptionTable bad}];

-1 "PASS test_calibration_validation";
