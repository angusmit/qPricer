/ test_calibration_bucket_report.q
\l core/init.q
spotTable:([] underlying:enlist `AAPL; spot:enlist 100f);
volTable:([] underlying:enlist `AAPL; volatility:enlist 0.2);
rateTable:([] expiry:enlist 1f; riskFreeRate:enlist 0.05);
divTable:([] underlying:enlist `AAPL; dividendYield:enlist 0f);
mktBook:.marketbook.createMarketDataBook[spotTable;volTable;rateTable;divTable];
p1:.validation.blackScholesClosedForm[`call;100f;100f;1f;0.05;0f;0.2];
optionTable:enlist `optionId`underlying`optionType`strike`expiry`marketPrice!(1;`AAPL;`call;100f;1f;p1);

bsCalib:.calibration.calibrateBlackScholesVolSurface[optionTable;mktBook;()!()];
compTable:.modelcompare.priceWithModel[`blackScholesSurface;optionTable;mktBook;bsCalib;()!()];

/ Maturity buckets
matBuckets:.calibreport.maturityBucketReport compTable;
.testutil.assertTrue[(count matBuckets)>0;"maturity buckets exist"];
.testutil.assertTrue[`maturityBucket in cols matBuckets;"has maturityBucket"];

/ Moneyness buckets
monBuckets:.calibreport.moneynessBucketReport[compTable;mktBook];
.testutil.assertTrue[(count monBuckets)>0;"moneyness buckets exist"];
.testutil.assertTrue[`moneynessBucket in cols monBuckets;"has moneynessBucket"];

/ Strike buckets
strikeBuckets:.calibreport.strikeBucketReport compTable;
.testutil.assertTrue[(count strikeBuckets)>0;"strike buckets exist"];

-1 "PASS test_calibration_bucket_report: matBuckets=",string[count matBuckets],", monBuckets=",string count monBuckets;
