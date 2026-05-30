\l core/init.q
pathCfg:`names`weights`spot0`vols`correlationMatrix`steps`stepYears`riskFreeRate`dividendYields`seed!(
    `A`B`C;1f%3 3 3f;100 80 120f;0.20 0.30 0.25;
    (1 0.4 0.2f;0.4 1 0.3f;0.2 0.3 1f);
    5;1f%252f;0.02;0 0 0f;7);
bundle:.strategy.path.fromCorrelated pathCfg;
requiredBundleKeys:`names`weights`pathTables`indexPath`correlationMatrix`vols;
.testutil.assertTrue[all requiredBundleKeys in key bundle;"bundle keys present"];
indexPath:bundle`indexPath;
schemaCols:`stepIndex`stepDate`spot`volatility`riskFreeRate`dividendYield`marketPrice`status;
.testutil.assertTableColumns[indexPath;schemaCols;"indexPath standard schema"];
.testutil.assertTrue[5=count indexPath;"5 rows in indexPath"];
pathTables:bundle`pathTables;
.testutil.assertTrue[(`A`B`C)~asc key pathTables;"per-name tables present"];
.testutil.assertTrue[all 5=count each value pathTables;"each per-name table has 5 rows"];
/ Index spot at step 0 = weighted spot0
expectedIndexSpot0:sum (pathCfg`weights)*pathCfg`spot0;
.testutil.assertTrue[1e-10>abs (first indexPath`spot)-expectedIndexSpot0;"index spot at step 0 = weighted basket"];
/ Per-name first spots are inputs
.testutil.assertTrue[100f=first (pathTables`A)`spot;"spot0 A preserved"];
.testutil.assertTrue[80f=first (pathTables`B)`spot;"spot0 B preserved"];
.testutil.assertTrue[120f=first (pathTables`C)`spot;"spot0 C preserved"];
/ Index spot at each step = weighted sum of per-name spots
nSpots:count indexPath;
constructedIdxSpots:sum (pathCfg`weights)*pathTables[;`spot];
.testutil.assertTrue[1e-10>max abs (indexPath`spot)-constructedIdxSpots;"index spot path = weighted sum at every step"];
-1 "PASS test_strategy_path_correlated_schema: 3 names x 5 steps, basket identity holds";
