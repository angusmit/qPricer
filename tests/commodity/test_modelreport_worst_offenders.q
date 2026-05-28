\l lib/init.q
/ worstOffenders returns the top-N positions per metric, sorted by absolute
/ metric value descending. Verify per-metric row count, ranking, and ordering.

mcCfg:.montecarlo.defaultMcConfig[];
mcCfg:@[mcCfg;`pathCount;:;3000];
mcCfg:@[mcCfg;`timeStepCount;:;15];

setupCall:`optionType`strikePrice`expiry`riskFreeRate`forwardPrice`spotPrice!(`call;75f;1f;0.05;75f;75f);
setupPut:@[setupCall;(`optionType;`strikePrice);:;(`put;80f)];
setupCallItm:@[setupCall;`strikePrice;:;65f];
inputsAll:`black76`schwartz`schwartz2`mrjump!(
    enlist[`volatility]!enlist 0.30;
    `x0`params!(log 75f;`meanReversionSpeed`longRunLogMean`volatility!(1.5;log 75f;0.30));
    `shortFactor0`longFactor0`params!(0f;log 75f;`meanReversionSpeed`shortVolatility`longVolatility`longDrift`correlation!(1.5;0.30;0.15;0f;0.25));
    `x0`params`mcConfig!(log 75f;`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility!(2f;log 75f;0.30;1f;0.05;0.20);mcCfg));
scenCfg:`forwardShiftPct`spotShiftPct`volShift!(0.05;0.05;0.05);

pos1:`tradeId`commodity`optionSetup`modelInputs`scenarioConfig`quantity`contractMultiplier!(`A;`wti;setupCall;inputsAll;scenCfg;1f;1000f);
pos2:@[pos1;(`tradeId;`commodity;`optionSetup;`quantity);:;(`B;`brent;setupPut;-2f)];
pos3:@[pos1;(`tradeId;`commodity;`optionSetup;`quantity);:;(`C;`gas;setupCallItm;3f)];

positions:(pos1;pos2;pos3);
greekCfg:.commodity.modelreport.defaultGreekConfig[];
disCfg:.commodity.modelreport.defaultDisagreementConfig[];

portfolioResult:.commodity.modelreport.runPortfolioRisk[positions;greekCfg;disCfg];
offenders:.commodity.modelreport.worstOffenders[portfolioResult;2];

requiredCols:`tradeId`commodity`metricName`metricValue`severity`offenderRank;
.testutil.assertTableColumns[offenders;requiredCols;"worstOffenders schema"];

metricNamesSeries:distinct offenders`metricName;
expectedMetrics:`priceRangePct`scenarioPnlRange`primarySensitivityRange`volatilityVegaRange`jumpIntensitySensitivity;
.testutil.assertTrue[all expectedMetrics in metricNamesSeries;"all expected metricNames present"];
.testutil.assertTrue[10=count offenders;"five metrics x topN=2 -> 10 rows"];

countsByMetric:exec count i by metricName from offenders;
.testutil.assertTrue[all 2=value countsByMetric;"exactly 2 rows per metric"];

ranksByMetric:exec offenderRank by metricName from offenders;
.testutil.assertTrue[all (`long$1 2)~/:value ranksByMetric;"ranks 1..2 per metric"];

priceRangePctRows:offenders where (offenders`metricName)=`priceRangePct;
priceMetricValues:abs priceRangePctRows`metricValue;
.testutil.assertTrue[priceMetricValues[0]>=priceMetricValues 1;"priceRangePct rows sorted descending by abs metricValue"];

scenRows:offenders where (offenders`metricName)=`scenarioPnlRange;
scenMetricValues:abs scenRows`metricValue;
.testutil.assertTrue[scenMetricValues[0]>=scenMetricValues 1;"scenarioPnlRange rows sorted descending by abs metricValue"];

severitySymbols:distinct offenders`severity;
.testutil.assertTrue[all severitySymbols in `OK`warning;"severity values are OK or warning"];

topN1Offenders:.commodity.modelreport.worstOffenders[portfolioResult;1];
.testutil.assertTrue[5=count topN1Offenders;"topN=1 yields five rows total"];

-1 "PASS test_modelreport_worst_offenders: rows=",string[count offenders],", topN1Rows=",string[count topN1Offenders];
