/ test_calibration_model_ranking.q
\l lib/init.q
spotTable:([] underlying:enlist `AAPL; spot:enlist 100f);
volTable:([] underlying:enlist `AAPL; volatility:enlist 0.2);
rateTable:([] expiry:enlist 1f; riskFreeRate:enlist 0.05);
divTable:([] underlying:enlist `AAPL; dividendYield:enlist 0f);
mktBook:.marketbook.createMarketDataBook[spotTable;volTable;rateTable;divTable];
p1:.validation.blackScholesClosedForm[`call;100f;100f;1f;0.05;0f;0.2];
optionTable:enlist `optionId`underlying`optionType`strike`expiry`marketPrice!(1;`AAPL;`call;100f;1f;p1);

bsCalib:.calibration.calibrateBlackScholesVolSurface[optionTable;mktBook;()!()];
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
paramGrid:`initialVarianceList`longRunVarianceList`meanReversionList`volOfVolList`correlationList!(
    enlist 0.04;enlist 0.04;enlist 2.0;enlist 0.0;enlist 0.0);
mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(5000;25;42;0b;0b;0.95);
hestonCalib:.calibration.calibrateHestonGrid[optionTable;mkt;paramGrid;enlist[`mcConfig]!enlist mcConfig];

calibDict:`blackScholesSurface`hestonGrid!(bsCalib;hestonCalib);
compTable:.modelcompare.compareModels[`blackScholesSurface`hestonGrid;optionTable;mktBook;calibDict;enlist[`mcConfig]!enlist mcConfig];

ranked:.modelcompare.rankModels compTable;
.testutil.assertTrue[2=count ranked;"2 ranked models"];

ranks:ranked`rank;
.testutil.assertTrue[1 in ranks;"rank 1 exists"];

bestRow:.modelcompare.bestModel compTable;
.testutil.assertTrue[bestRow[`rank]=1;"best is rank 1"];
.testutil.assertTrue[not null bestRow`rmse;"best rmse finite"];

-1 "PASS test_calibration_model_ranking: best=",string[bestRow`modelName],", rmse=",string bestRow`rmse;
