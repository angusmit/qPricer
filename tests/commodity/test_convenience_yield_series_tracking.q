\l core/init.q
/ Part A: a synthetic curveHistory built from a KNOWN time-varying convenience
/ yield (varying X0 across dates at fixed kappa/vols drives the front slope).
/ Assert the recovered cy series tracks the known cy computed INDEPENDENTLY from
/ the generated market front slope.
kappa:1.0; shortVol:0.30; longVol:0.15; corrVal:0.30; rate:0.02; longFactor0:4.10;
tenors:0.08 0.25 0.5 0.75 1.0 1.5f;
dates:2020.01.06 2020.02.06 2020.03.06 2020.05.06 2020.08.06;
x0s:0.25 0.10 -0.05 -0.20 -0.30;
genOne:{[d;x0;tenors;kappaP;y0;shortVolP;longVolP;corrP]
    params:`meanReversionSpeed`shortVolatility`longVolatility`longDrift`correlation!(kappaP;shortVolP;longVolP;-0.05;corrP);
    prices:.commodity.schwartz2.futuresCurve[x0;y0;params;tenors];
    ([] asofDate:(count tenors)#d; tenor:tenors; price:prices; contractYM:(count tenors)#0; expiry:(count tenors)#d) };
hist:raze genOne[;;tenors;kappa;longFactor0;shortVol;longVol;corrVal]'[dates;x0s];

/ Known cy per date from the market front slope (independent of the fitter).
knownCy:{[d;histTbl;r] sub:`tenor xasc select tenor,price from histTbl where asofDate=d; p:sub`price; t:sub`tenor; r-((log p 1)-log p 0)%t[1]-t 0}[;hist;rate] each dates;

calCfg:`kappa`shortVolatility`longVolatility`correlation`riskFreeRate!(kappa;shortVol;longVol;corrVal;rate);
out:.commodity.curveCal.convenienceYieldSeries[hist;calCfg];
series:out`series;

.testutil.assertTrue[5=count series;"5 dates fitted"];
.testutil.assertTrue[(`asofDate`frontPrice`netConvenienceYield`slope`regime`fitRmse`X0`Y0`muY`status)~cols series;"series schema"];
.testutil.assertTrue[1e-6>max abs knownCy-series`netConvenienceYield;"recovered cy tracks known cy"];
.testutil.assertTrue[`backwardation=first series`regime;"early dates backwardation (cy > rate)"];
.testutil.assertTrue[`contango=last series`regime;"late dates contango (cy < rate)"];
.testutil.assertTrue[1=count out`transitions;"exactly one regime transition"];
tr:first out`transitions;
.testutil.assertTrue[`backwardation`contango~tr`fromRegime`toRegime;"transition backwardation -> contango"];
.testutil.assertTrue[0=count out`skipped;"no skipped dates (all prices positive)"];

-1 "PASS test_convenience_yield_series_tracking: maxErr=",string max abs knownCy-series`netConvenienceYield;
