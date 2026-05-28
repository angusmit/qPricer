/ test_portfolio_lookback_products.q
\l lib/init.q

/ Test lookback via direct pricing (portfolio routing needs flat market data)
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
mcConfig:.montecarlo.defaultMcConfig[];

/ Fixed call lookback
fixedCall:`tradeId`underlying`productType`exerciseStyle`optionType`lookbackStyle`strike`expiry`notional`observationCount!(
    1;`AAPL;`lookbackOption;`european;`call;`fixed;100f;1f;100000f;50);
fixedResult:.lookback.priceLookbackOption[fixedCall;mkt;enlist[`mcConfig]!enlist mcConfig];
.testutil.assertTrue[fixedResult[`status]~`OK;"fixed call OK"];
.testutil.assertTrue[fixedResult[`unitPrice]>0f;"fixed call price positive"];
.testutil.assertNear[fixedResult`notionalPrice;fixedResult[`unitPrice]*100000f;1f;"notional correct"];

/ American lookback should error
americanLB:`tradeId`underlying`productType`exerciseStyle`optionType`lookbackStyle`strike`expiry`notional`observationCount!(
    2;`AAPL;`lookbackOption;`american;`call;`fixed;100f;1f;100000f;50);
americanErr:@[{.lookback.priceLookbackOption[x;mkt;enlist[`mcConfig]!enlist mcConfig];`NO_ERROR};americanLB;{`ERROR}];
.testutil.assertTrue[americanErr~`ERROR;"American lookback errors"];

/ Floating put
floatingPut:`tradeId`underlying`productType`exerciseStyle`optionType`lookbackStyle`strike`expiry`notional`observationCount!(
    3;`AAPL;`lookbackOption;`european;`put;`floating;0Nf;1f;100000f;50);
floatResult:.lookback.priceLookbackOption[floatingPut;mkt;enlist[`mcConfig]!enlist mcConfig];
.testutil.assertTrue[floatResult[`status]~`OK;"floating put OK"];
.testutil.assertTrue[floatResult[`unitPrice]>0f;"floating put price positive"];

-1 "PASS test_portfolio_lookback_products: fixedCall=",string[fixedResult`unitPrice],", floatingPut=",string floatResult`unitPrice;
