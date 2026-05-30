\l core/init.q
.test.expectError:{[testName;fn]
    testResult:@[{x[];`NO_ERROR};fn;{`ERROR}];
    if[testResult~`NO_ERROR; '"Expected error but got success: ",testName];
    -1 "  PASS ",testName;
 };

.test.expectError["empty PnL";{.riskdist.validatePnlVector `float$()}];
.test.expectError["null PnL";{.riskdist.validatePnlVector 1 0N 3f}];
.test.expectError["empty scenario";{.riskdist.validateScenarioResult ()}];
.test.expectError["missing scenarioName";{.riskdist.validateScenarioResult enlist `pnl!(enlist 1f)}];
.test.expectError["confidence 0";{.var.validateConfidenceLevel 0f}];
.test.expectError["confidence 1";{.var.validateConfidenceLevel 1f}];
.test.expectError["confidence negative";{.var.validateConfidenceLevel -0.5}];

-1 "PASS test_var_validation";
