\l lib/init.q
/ validateOptionSetup rejects malformed setups; per-model pricing wrappers
/ degrade to an ERROR row when their inputs are missing instead of crashing
/ the whole report; comparisonSummary aggregates only OK rows.

baseSetup:`optionType`strikePrice`expiry`riskFreeRate`forwardPrice`spotPrice!(`call;75f;1f;0.05;75f;75f);

badType:@[baseSetup;`optionType;:;`straddle];
r1:@[.commodity.modelreport.validateOptionSetup;badType;{`ERROR}];
.testutil.assertTrue[r1~`ERROR;"unsupported optionType rejected"];

badStrike:@[baseSetup;`strikePrice;:;-10f];
r2:@[.commodity.modelreport.validateOptionSetup;badStrike;{`ERROR}];
.testutil.assertTrue[r2~`ERROR;"negative strike rejected"];

badForward:@[baseSetup;`forwardPrice;:;0f];
r3:@[.commodity.modelreport.validateOptionSetup;badForward;{`ERROR}];
.testutil.assertTrue[r3~`ERROR;"zero forwardPrice rejected"];

missingKey:`optionType`strikePrice`expiry`riskFreeRate`forwardPrice!(`call;75f;1f;0.05;75f);
r4:@[.commodity.modelreport.validateOptionSetup;missingKey;{`ERROR}];
.testutil.assertTrue[r4~`ERROR;"missing spotPrice rejected"];

/ Missing model inputs cause per-row failure, not full crash.
partialInputs:`black76`schwartz!(
    enlist[`volatility]!enlist 0.30;
    `x0`params!(log 75f;`meanReversionSpeed`longRunLogMean`volatility!(1.5;log 75f;0.30)));

partialPriceTable:.commodity.modelreport.priceAllModels[baseSetup;partialInputs];
.testutil.assertTrue[4=count partialPriceTable;"partial inputs still produce four rows"];

statusByModel:exec pricingStatus by modelName from partialPriceTable;
.testutil.assertTrue[`OK=first statusByModel`black76;"black76 OK with valid inputs"];
.testutil.assertTrue[`OK=first statusByModel`schwartz;"schwartz OK with valid inputs"];
.testutil.assertTrue[`ERROR=first statusByModel`schwartz2;"schwartz2 ERROR row for missing input"];
.testutil.assertTrue[`ERROR=first statusByModel`mrjump;"mrjump ERROR row for missing input"];

partialSummary:.commodity.modelreport.comparisonSummary partialPriceTable;
.testutil.assertTrue[4=partialSummary`modelCount;"modelCount counts all rows"];
.testutil.assertTrue[2=partialSummary`okModelCount;"okModelCount counts OK rows only"];
.testutil.assertTrue[`OK=partialSummary`status;"summary status OK while at least one model priced"];

/ All-fail case must still produce a usable summary with status ERROR.
emptyInputs:()!();
allFailTable:.commodity.modelreport.priceAllModels[baseSetup;emptyInputs];
.testutil.assertTrue[all (allFailTable`pricingStatus)=`ERROR;"all model rows ERROR when no inputs supplied"];
allFailSummary:.commodity.modelreport.comparisonSummary allFailTable;
.testutil.assertTrue[`ERROR=allFailSummary`status;"summary status ERROR when no model priced"];
.testutil.assertTrue[0=allFailSummary`okModelCount;"okModelCount = 0 in all-fail case"];

-1 "PASS test_modelreport_validation: partialOk=",string[partialSummary`okModelCount],"/",string partialSummary`modelCount;
