\l core/init.q
.test.expectError:{[testName;fn]
    testResult:@[{x[];`NO_ERROR};fn;{`ERROR}];
    if[testResult~`NO_ERROR; '"Expected error: ",testName];
    -1 "  PASS ",testName;
 };

.test.expectError["empty limit table";{.limits.validateLimitTable ()}];
.test.expectError["empty metric table";{
    lt:enlist `limitId`scopeType`scopeValue`metricName`limitValue`warningPct`hardLimitPct`direction`enabled!(1;`portfolio;`ALL;`VaR95;100f;0.8;1.0;`lessThan;1b);
    .limits.evaluateLimits[lt;()]}];

/ Missing metric returns error severity
metricTable:enlist `scopeType`scopeValue`metricName`metricValue`source!(`portfolio;`ALL;`VaR95;100f;`var);
lt:enlist `limitId`scopeType`scopeValue`metricName`limitValue`warningPct`hardLimitPct`direction`enabled!(1;`portfolio;`ALL;`UnknownMetric;100f;0.8;1.0;`lessThan;1b);
checkResult:.limits.evaluateLimits[lt;metricTable];
.testutil.assertTrue[(checkResult 0)[`severity]=`error;"unknown metric -> error"];

-1 "PASS test_limit_monitoring_validation";
