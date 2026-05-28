\l lib/init.q
/ greeksBlack76 produces FD delta and vega per unit-of-bump. Compare against
/ closed-form Greeks from .commodity.black76.greeks: closed-form vega is per
/ 1% vol move (divided by 100), so multiplying it by 100 recovers the
/ per-unit-volatility derivative that the FD path produces.

optionSetup:`optionType`strikePrice`expiry`riskFreeRate`forwardPrice`spotPrice!(`call;75f;1f;0.05;75f;75f);
modelInputs:enlist[`black76]!enlist enlist[`volatility]!enlist 0.30;

greekCfg:.commodity.modelreport.defaultGreekConfig[];
greekTbl:.commodity.modelreport.greeksBlack76[optionSetup;modelInputs;greekCfg];

requiredCols:`modelName`basePrice`forwardDelta`volatilityVega`pricingStatus`pricingErrorMessage;
.testutil.assertTableColumns[greekTbl;requiredCols;"greeksBlack76 schema"];
.testutil.assertTrue[1=count greekTbl;"greeksBlack76 returns one row"];

greekRow:greekTbl 0;
.testutil.assertTrue[greekRow[`pricingStatus]=`OK;"black76 greeks status OK"];
.testutil.assertTrue[greekRow[`basePrice]>0f;"basePrice positive for ITM-region call"];
.testutil.assertTrue[not null greekRow`forwardDelta;"forwardDelta finite"];
.testutil.assertTrue[greekRow[`forwardDelta]>0f;"forwardDelta positive for call"];
.testutil.assertTrue[greekRow[`forwardDelta]<1f;"forwardDelta below 1 for non-deep call"];
.testutil.assertTrue[greekRow[`volatilityVega]>0f;"volatilityVega positive for call"];

closedForm:.commodity.black76.greeks[`call;75f;75f;1f;0.30;0.05];
closedDeltaVal:closedForm`delta;
closedVegaRaw:100f*closedForm`vega;
.testutil.assertNear[greekRow`forwardDelta;closedDeltaVal;1e-4;"FD forwardDelta within 1e-4 of closed-form delta"];
.testutil.assertNear[greekRow`volatilityVega;closedVegaRaw;1e-2;"FD vega within 1e-2 of closed-form vega-per-unit"];

-1 "PASS test_modelreport_greeks_black76: basePrice=",string[greekRow`basePrice],", delta=",string[greekRow`forwardDelta]," (closed ",string[closedDeltaVal],"), vega=",string[greekRow`volatilityVega]," (closed ",string[closedVegaRaw],")";
