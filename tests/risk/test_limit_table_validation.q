\l core/init.q
/ Valid table
validTable:();
validTable:validTable,enlist `limitId`scopeType`scopeValue`metricName`limitValue`warningPct`hardLimitPct`direction`enabled!(1;`portfolio;`ALL;`VaR95;1000000f;0.8;1.0;`lessThan;1b);
validTable:validTable,enlist `limitId`scopeType`scopeValue`metricName`limitValue`warningPct`hardLimitPct`direction`enabled!(2;`book;`EQD;`DeltaCash;500000f;0.8;1.0;`absLessThan;1b);
.limits.validateLimitTable validTable;
-1 "  PASS valid table";

/ Duplicate limitId
dupTable:();
dupTable:dupTable,enlist `limitId`scopeType`scopeValue`metricName`limitValue`warningPct`hardLimitPct`direction`enabled!(1;`portfolio;`ALL;`VaR95;1000f;0.8;1.0;`lessThan;1b);
dupTable:dupTable,enlist `limitId`scopeType`scopeValue`metricName`limitValue`warningPct`hardLimitPct`direction`enabled!(1;`book;`EQD;`Delta;500f;0.8;1.0;`lessThan;1b);
r:@[.limits.validateLimitTable;dupTable;{`ERROR}];
.testutil.assertTrue[r~`ERROR;"duplicate limitId fails"];

/ Negative limitValue
negTable:enlist `limitId`scopeType`scopeValue`metricName`limitValue`warningPct`hardLimitPct`direction`enabled!(1;`portfolio;`ALL;`VaR95;-100f;0.8;1.0;`lessThan;1b);
r2:@[.limits.validateLimitTable;negTable;{`ERROR}];
.testutil.assertTrue[r2~`ERROR;"negative limitValue fails"];

/ warningPct > hardLimitPct
warnTable:enlist `limitId`scopeType`scopeValue`metricName`limitValue`warningPct`hardLimitPct`direction`enabled!(1;`portfolio;`ALL;`VaR95;100f;1.2;1.0;`lessThan;1b);
r3:@[.limits.validateLimitTable;warnTable;{`ERROR}];
.testutil.assertTrue[r3~`ERROR;"warningPct > hardLimitPct fails"];

/ Invalid direction
badDir:enlist `limitId`scopeType`scopeValue`metricName`limitValue`warningPct`hardLimitPct`direction`enabled!(1;`portfolio;`ALL;`VaR95;100f;0.8;1.0;`badDir;1b);
r4:@[.limits.validateLimitTable;badDir;{`ERROR}];
.testutil.assertTrue[r4~`ERROR;"invalid direction fails"];

-1 "PASS test_limit_table_validation";
