/ test_model_comparison.q
\l lib/init.q
spotTable:([] underlying:enlist `AAPL; spot:enlist 100f);
volTable:([] underlying:enlist `AAPL; volatility:enlist 0.2);
rateTable:([] expiry:enlist 1f; riskFreeRate:enlist 0.05);
divTable:([] underlying:enlist `AAPL; dividendYield:enlist 0f);
mktBook:.marketbook.createMarketDataBook[spotTable;volTable;rateTable;divTable];

p1:.validation.blackScholesClosedForm[`call;100f;100f;1f;0.05;0f;0.2];
p2:.validation.blackScholesClosedForm[`put;100f;100f;1f;0.05;0f;0.2];
optionTable:();
optionTable:optionTable,enlist `optionId`underlying`optionType`strike`expiry`marketPrice!(1;`AAPL;`call;100f;1f;p1);
optionTable:optionTable,enlist `optionId`underlying`optionType`strike`expiry`marketPrice!(2;`AAPL;`put;100f;1f;p2);

/ BS calibration
bsCalib:.calibration.calibrateBlackScholesVolSurface[optionTable;mktBook;()!()];

/ Small Heston grid calibration
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
paramGrid:`initialVarianceList`longRunVarianceList`meanReversionList`volOfVolList`correlationList!(
    enlist 0.04;enlist 0.04;enlist 2.0;0.0 0.3;enlist -0.7);
mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(5000;25;42;0b;0b;0.95);
hestonCalib:.calibration.calibrateHestonGrid[optionTable;mkt;paramGrid;enlist[`mcConfig]!enlist mcConfig];

calibDict:`blackScholesSurface`hestonGrid!(bsCalib;hestonCalib);

compTable:.modelcompare.compareModels[`blackScholesSurface`hestonGrid;optionTable;mktBook;calibDict;enlist[`mcConfig]!enlist mcConfig];

.testutil.assertTrue[(count compTable)>0;"comparison has rows"];
.testutil.assertTrue[`modelName in cols compTable;"has modelName column"];
.testutil.assertTrue[`modelPrice in cols compTable;"has modelPrice column"];

okCount:sum (compTable`status)=`OK;
.testutil.assertTrue[okCount>0;"some rows priced OK"];

-1 "PASS test_model_comparison: rows=",string[count compTable],", okRows=",string okCount;
