\l core/init.q
.test.expectError:{[testName;fn]
    testResult:@[{x[];`NO_ERROR};fn;{`ERROR}];
    if[testResult~`NO_ERROR; '"Expected error but got success: ",testName];
    -1 "  PASS ",testName;
 };
.test.expectError["empty shock table";{.histscen.validateHistoricalShockTable ()}];
.test.expectError["missing column";{.histscen.validateHistoricalShockTable enlist `scenarioDate`eventName!(2020.01.01;`test)}];
.test.expectError["spotReturn <= -1";{
    bad:enlist `scenarioDate`eventName`underlying`spotReturn`volatilityShift`rateShift`dividendShift!(2020.01.01;`test;`AAPL;-1.5;0f;0f;0f);
    .histscen.validateHistoricalShockTable bad}];
-1 "PASS test_historical_replay_validation";
