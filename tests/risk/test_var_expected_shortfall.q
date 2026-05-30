\l core/init.q
pnlVec:-10 -5 0 2 4 6f;

varVal:.var.valueAtRisk[pnlVec;0.95];
.testutil.assertTrue[varVal>0f;"VaR is positive loss"];

esVal:.var.expectedShortfall[pnlVec;0.95];
.testutil.assertTrue[esVal>=varVal;"ES >= VaR"];

/ Multiple confidence levels
var95:.var.valueAtRisk[pnlVec;0.95];
var90:.var.valueAtRisk[pnlVec;0.9];
.testutil.assertTrue[var95>=var90;"higher confidence -> higher VaR"];

/ VaR/ES combined
vesResult:.var.varExpectedShortfall[pnlVec;0.95];
.testutil.assertTrue[vesResult[`status]~`OK;"status OK"];
.testutil.assertTrue[vesResult[`observationCount]=6;"6 observations"];

/ Invalid confidence
invalidResult:@[.var.validateConfidenceLevel;0f;{`ERROR}];
.testutil.assertTrue[invalidResult~`ERROR;"0 confidence fails"];
invalidResult2:@[.var.validateConfidenceLevel;1f;{`ERROR}];
.testutil.assertTrue[invalidResult2~`ERROR;"1 confidence fails"];

/ Empty vector
emptyResult:@[.var.valueAtRisk[;0.95];`float$();{`ERROR}];
.testutil.assertTrue[emptyResult~`ERROR;"empty vector fails"];

-1 "PASS test_var_expected_shortfall: VaR95=",string[varVal],", ES95=",string esVal;
