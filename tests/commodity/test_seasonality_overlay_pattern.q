\l core/init.q
/ A seasonal curve shows winter > summer at the right tenors. With phaseYears=0
/ the peak is at integer-year delivery and the trough at the half-year. At step 0
/ the delivery time equals the tenor and (contango=0) the baseline is flat, so the
/ seasonal multiplier exp(factor) is the only term shaping the curve.
seasonCfg:`amplitude`phaseYears!(0.2;0f);
base:`tenors`evolutionModel`evolutionParams`steps`stepYears`seed`seasonCfg!(
    0.04 0.5 0.96f;`simple;
    `spot0`drift`volatility`contango`riskFreeRate!(60f;0f;0.2;0f;0.02);
    6;1f%252f;11;seasonCfg);
bundle:.strategy.path.fromFuturesCurve base;
snaps:bundle`curveSnapshots;
winterNear:first exec futuresPrice from snaps where stepIndex=0, tenor=0.04;
summerMid: first exec futuresPrice from snaps where stepIndex=0, tenor=0.5;
winterFar: first exec futuresPrice from snaps where stepIndex=0, tenor=0.96;
.testutil.assertTrue[winterNear>summerMid;"winter (t~0) forward > summer (t=0.5)"];
.testutil.assertTrue[winterFar>summerMid;"winter (t~1) forward > summer (t=0.5)"];

/ Known-answer: the winter/summer ratio equals exp(factor diff) (flat baseline).
factorWinter:.commodity.seasonality.factor[0.04;seasonCfg];
factorSummer:.commodity.seasonality.factor[0.5;seasonCfg];
expectedRatio:exp factorWinter-factorSummer;
.testutil.assertNear[winterNear%summerMid;expectedRatio;1e-9;"seasonal ratio = exp(factor difference)"];

-1 "PASS test_seasonality_overlay_pattern: winterNear=",(string winterNear),", summerMid=",string summerMid;
