\l lib/init.q
/ runComparisonWithGreeks is the top-level convenience wrapper: it returns the
/ runComparison report (basePrices, priceDifferences, summary, scenarioPrices,
/ scenarioPnL) joined with greeksTable and greeksSummary. Verify the resulting
/ dictionary exposes all seven keys and that the embedded tables are non-empty.

optionSetup:`optionType`strikePrice`expiry`riskFreeRate`forwardPrice`spotPrice!(`call;75f;1f;0.05;75f;75f);

mcCfg:.montecarlo.defaultMcConfig[];
mcCfg:@[mcCfg;`pathCount;:;8000];
mcCfg:@[mcCfg;`timeStepCount;:;25];

modelInputs:`black76`schwartz`schwartz2`mrjump!(
    enlist[`volatility]!enlist 0.30;
    `x0`params!(log 75f;`meanReversionSpeed`longRunLogMean`volatility!(1.5;log 75f;0.30));
    `shortFactor0`longFactor0`params!(0f;log 75f;`meanReversionSpeed`shortVolatility`longVolatility`longDrift`correlation!(1.5;0.30;0.15;0f;0.25));
    `x0`params`mcConfig!(log 75f;`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility!(2f;log 75f;0.30;1f;0.05;0.20);mcCfg));

scenarioCfg:`forwardShiftPct`spotShiftPct`volShift!(0.05;0.05;0.05);
greekCfg:.commodity.modelreport.defaultGreekConfig[];

combinedReport:.commodity.modelreport.runComparisonWithGreeks[optionSetup;modelInputs;scenarioCfg;1f;1000f;greekCfg];

expectedKeys:`basePrices`priceDifferences`summary`scenarioPrices`scenarioPnL`greeksTable`greeksSummary;
.testutil.assertTrue[all expectedKeys in key combinedReport;"combinedReport has all expected keys"];

basePricesTbl:combinedReport`basePrices;
.testutil.assertTrue[4=count basePricesTbl;"basePrices has four rows"];

greeksTbl:combinedReport`greeksTable;
.testutil.assertTrue[4=count greeksTbl;"greeksTable has four rows"];

greekSummary:combinedReport`greeksSummary;
.testutil.assertTrue[`OK=greekSummary`status;"greeksSummary status OK"];
.testutil.assertTrue[4=greekSummary`okModelCount;"all four model greeks OK"];

pnLTbl:combinedReport`scenarioPnL;
.testutil.assertTrue[4=count pnLTbl;"scenarioPnL still has four rows"];
.testutil.assertTrue[all (pnLTbl`status)=`OK;"all scenarioPnL rows OK"];

-1 "PASS test_modelreport_comparison_with_greeks: greekOkCount=",string[greekSummary`okModelCount]," maxPrimaryDelta=",string[greekSummary`maxPrimaryDelta];
