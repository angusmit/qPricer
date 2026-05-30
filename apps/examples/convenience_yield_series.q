\l core/init.q
/ ============================================================================
/ convenience_yield_series.q - real-data example (NOT a test).
/ ----------------------------------------------------------------------------
/ Part A across the REAL WTI curve history 2020-2021: calibrate the curve per
/ as-of date (kappa fixed; a single snapshot cannot identify it) to produce a
/ convenience-yield series, telling the backwardation -> COVID super-contango ->
/ recovery story. Dates with non-positive prices (the April-2020 negative-settle
/ window) are auto-skipped. CSVs are user-supplied and gitignored.
/ ============================================================================
crudeDir:"data/barchart/CRUDE";
longTable:.parser.crude.loadAll crudeDir;
allDates:asc distinct longTable`date;
inRange:allDates where (allDates>=2020.01.01) and allDates<=2021.12.31;
asofDates:inRange where 0=(til count inRange) mod 21;
-1 "Building curve history over ",(string count asofDates)," monthly-ish as-of dates ...";
hist:.parser.crude.curveHistory[longTable;asofDates];

calCfg:`kappa`shortVolatility`longVolatility`correlation`riskFreeRate!(1.0;0.35;0.20;0.30;0.02);
out:.commodity.curveCal.convenienceYieldSeries[hist;calCfg];
series:out`series;

-1 "";
-1 "Convenience-yield series (kappa fixed at 1.0; Part B identifies it from dynamics):";
show select asofDate,frontPrice,netConvenienceYield,regime,fitRmse from series;
-1 "";
-1 "Regime transitions (backwardation <-> contango flips):";
show out`transitions;
-1 "";
-1 "Skipped dates (non-positive prices, incl. April-2020 negative settles): ",-3!out`skipped;
-1 "Story: WTI sat in backwardation pre-COVID (convenience yield > carry), flipped";
-1 "to steep contango in the 2020 demand collapse, then recovered toward backwardation.";
