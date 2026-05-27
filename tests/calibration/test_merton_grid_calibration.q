/ test_merton_grid_calibration.q
\l lib/init.q
/ Synthetic prices from known Merton params
trueParams:`volatility`jumpIntensity`jumpMean`jumpVolatility`riskFreeRate`dividendYield!(0.2;0.5;-0.1;0.2;0.05;0.0);
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);

p1:.merton.priceEuropeanSeries[`call;100f;90f;1f;trueParams;30];
p2:.merton.priceEuropeanSeries[`call;100f;100f;1f;trueParams;30];
p3:.merton.priceEuropeanSeries[`put;100f;110f;1f;trueParams;30];

optionTable:();
optionTable:optionTable,enlist `optionId`underlying`optionType`strike`expiry`marketPrice!(1;`AAPL;`call;90f;1f;p1);
optionTable:optionTable,enlist `optionId`underlying`optionType`strike`expiry`marketPrice!(2;`AAPL;`call;100f;1f;p2);
optionTable:optionTable,enlist `optionId`underlying`optionType`strike`expiry`marketPrice!(3;`AAPL;`put;110f;1f;p3);

paramGrid:`volatilityList`jumpIntensityList`jumpMeanList`jumpVolatilityList!(
    0.15 0.2 0.25;
    0.0 0.5 1.0;
    -0.1 0.0;
    0.1 0.2 0.3);

calibTable:.merton.calibrateGrid[optionTable;mkt;paramGrid;()!()];
.testutil.assertTrue[(count calibTable)>0;"calibration has rows"];

bestRow:.merton.bestCalibration calibTable;
.testutil.assertTrue[not null bestRow`rmse;"best rmse finite"];
.testutil.assertTrue[bestRow[`pricedRows]>0;"best has priced rows"];
.testutil.assertTrue[bestRow[`rmse]<0.1;"best rmse reasonable"];

-1 "PASS test_merton_grid_calibration: gridRows=",string[count calibTable],", bestRmse=",string[bestRow`rmse],", bestVol=",string bestRow`volatility;
