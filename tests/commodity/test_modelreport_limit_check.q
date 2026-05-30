\l core/init.q
/ checkLimit returns dict with status OK / warning / breach / ERROR depending
/ on metric value vs warning/breach thresholds. Verify all five branches:
/ OK, warning, breach (with breachAmount), null metric (ERROR), invalid
/ thresholds where warning > breach (ERROR).

okResult:.commodity.modelreport.checkLimit[`grossPriceRangeExposure;100f;5000f;10000f];
.testutil.assertTrue[`OK=okResult`status;"OK case status OK"];
.testutil.assertTrue[0f=okResult`breachAmount;"OK case breachAmount = 0"];

warnResult:.commodity.modelreport.checkLimit[`grossPriceRangeExposure;7500f;5000f;10000f];
.testutil.assertTrue[`warning=warnResult`status;"warning case status warning"];
.testutil.assertTrue[0f=warnResult`breachAmount;"warning case breachAmount = 0"];

breachResult:.commodity.modelreport.checkLimit[`grossPriceRangeExposure;15000f;5000f;10000f];
.testutil.assertTrue[`breach=breachResult`status;"breach case status breach"];
.testutil.assertNear[breachResult`breachAmount;5000f;1e-12;"breach case breachAmount = metric - breach"];
.testutil.assertTrue[0<count breachResult`message;"breach case message non-empty"];

nullResult:.commodity.modelreport.checkLimit[`grossPriceRangeExposure;0Nf;5000f;10000f];
.testutil.assertTrue[`ERROR=nullResult`status;"null metric case status ERROR"];
.testutil.assertTrue["metricValue is null"~nullResult`message;"null metric error message specific"];

invalidResult:.commodity.modelreport.checkLimit[`grossPriceRangeExposure;100f;10000f;5000f];
.testutil.assertTrue[`ERROR=invalidResult`status;"invalid threshold (warning > breach) status ERROR"];
.testutil.assertTrue["Invalid thresholds"~invalidResult`message;"invalid threshold error message specific"];

nullThr:.commodity.modelreport.checkLimit[`grossPriceRangeExposure;100f;0Nf;10000f];
.testutil.assertTrue[`ERROR=nullThr`status;"null warningThreshold rejected"];

equalThr:.commodity.modelreport.checkLimit[`grossPriceRangeExposure;5000f;5000f;5000f];
.testutil.assertTrue[`breach=equalThr`status;"metric equal to breach threshold is breach"];

-1 "PASS test_modelreport_limit_check: breachAmt=",string[breachResult`breachAmount];
