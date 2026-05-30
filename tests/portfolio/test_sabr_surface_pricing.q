/ test_sabr_surface_pricing.q
\l core/init.q

mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
sabrCalib:`alpha`beta`rho`nu!(0.2;0.5;0.0;0.0001);

/ European call
callOption:`optionId`underlying`optionType`strike`expiry`notional!(1;`AAPL;`call;100f;1f;100000f);
callResult:.sabr.priceEuropeanWithSabr[callOption;mkt;sabrCalib;()!()];
.testutil.assertTrue[callResult[`status]~`OK;"call pricing OK"];
.testutil.assertTrue[callResult[`unitPrice]>0f;"call price positive"];
.testutil.assertTrue[callResult[`impliedVolatility]>0f;"IV positive"];

/ European put
putOption:`optionId`underlying`optionType`strike`expiry`notional!(2;`AAPL;`put;100f;1f;100000f);
putResult:.sabr.priceEuropeanWithSabr[putOption;mkt;sabrCalib;()!()];
.testutil.assertTrue[putResult[`status]~`OK;"put pricing OK"];
.testutil.assertTrue[putResult[`unitPrice]>0f;"put price positive"];

-1 "PASS test_sabr_surface_pricing: callPrice=",string[callResult`unitPrice],", putPrice=",string putResult`unitPrice;
