/ test_calibration_report_export.q
\l lib/init.q
spotTable:([] underlying:enlist `AAPL; spot:enlist 100f);
volTable:([] underlying:enlist `AAPL; volatility:enlist 0.2);
rateTable:([] expiry:enlist 1f; riskFreeRate:enlist 0.05);
divTable:([] underlying:enlist `AAPL; dividendYield:enlist 0f);
mktBook:.marketbook.createMarketDataBook[spotTable;volTable;rateTable;divTable];
p1:.validation.blackScholesClosedForm[`call;100f;100f;1f;0.05;0f;0.2];
optionTable:enlist `optionId`underlying`optionType`strike`expiry`marketPrice!(1;`AAPL;`call;100f;1f;p1);

bsCalib:.calibration.calibrateBlackScholesVolSurface[optionTable;mktBook;()!()];
compTable:.modelcompare.priceWithModel[`blackScholesSurface;optionTable;mktBook;bsCalib;()!()];

exportResult:.calibreport.exportCalibrationReports[compTable;".";`testCalib];
.testutil.assertTrue[(count exportResult)>0;"export results exist"];
exportRow:exportResult 0;
.testutil.assertTrue[`reportName in key exportRow;"has reportName"];
.testutil.assertTrue[`status in key exportRow;"has status"];

-1 "PASS test_calibration_report_export: reports=",string count exportResult;
