\l core/init.q
/ Part A: dates whose curve contains a non-positive price are SKIPPED (excludes
/ the April-2020 negative-settle window), the rest are fitted, and regimes are
/ classified correctly. No crash.
kappa:1.0; shortVol:0.30; longVol:0.15; corrVal:0.30; rate:0.02; longFactor0:4.10;
tenors:0.08 0.25 0.5 0.75 1.0 1.5f;
params:`meanReversionSpeed`shortVolatility`longVolatility`longDrift`correlation!(kappa;shortVol;longVol;-0.05;corrVal);
backwardatedPrices:.commodity.schwartz2.futuresCurve[0.20;longFactor0;params;tenors];
contangoPrices:.commodity.schwartz2.futuresCurve[-0.20;longFactor0;params;tenors];
mk:{[d;prices;tenors] ([] asofDate:(count tenors)#d; tenor:tenors; price:prices; contractYM:(count tenors)#0; expiry:(count tenors)#d)};

/ Middle date carries a negative price (April-2020 style).
negPrices:backwardatedPrices;
negPrices[2]:-5.0;
hist:raze (mk[2020.01.06;backwardatedPrices;tenors];mk[2020.04.20;negPrices;tenors];mk[2020.08.06;contangoPrices;tenors]);

calCfg:`kappa`shortVolatility`longVolatility`correlation`riskFreeRate!(kappa;shortVol;longVol;corrVal;rate);
out:.commodity.curveCal.convenienceYieldSeries[hist;calCfg];
series:out`series;

.testutil.assertTrue[2=count series;"2 dates fitted (negative-price date skipped)"];
.testutil.assertTrue[2020.04.20 in out`skipped;"April negative-price date recorded as skipped"];
.testutil.assertTrue[not 2020.04.20 in series`asofDate;"skipped date absent from the series"];
.testutil.assertTrue[`backwardation=first series`regime;"first date classified backwardation (X0>0)"];
.testutil.assertTrue[`contango=last series`regime;"last date classified contango (X0<0)"];

-1 "PASS test_convenience_yield_series_skip: skipped=",-3!out`skipped;
