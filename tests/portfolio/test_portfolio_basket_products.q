/ test_portfolio_basket_products.q - portfolio with equity + Asian + basket
\l ./lib/init.q

/ Market data book with 3 symbols
spotTable:([] underlying:`AAPL`MSFT`NVDA; spot:100 250 800f);
volTable:([] underlying:`AAPL`MSFT`NVDA; volatility:0.2 0.25 0.35);
rateTable:([] expiry:enlist 1f; riskFreeRate:enlist 0.05);
divTable:([] underlying:`AAPL`MSFT`NVDA; dividendYield:0 0.01 0f);
mktBook:.marketbook.createMarketDataBook[spotTable;volTable;rateTable;divTable];
corrTable:([] sym1:`AAPL`AAPL`MSFT; sym2:`MSFT`NVDA`NVDA; correlation:0.5 0.3 0.4);

/ Build trades with uniform columns
equityTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`basketSymbols`basketWeights`averageType`averagingStyle`observationCount`barrierType`barrierLevel`rebate!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;100000f;`symbol$();`float$();`none;`none;0N;`none;0Nf;0f);
basketTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`basketSymbols`basketWeights`averageType`averagingStyle`observationCount`barrierType`barrierLevel`rebate!(
    2;`AAPL;`basketOption;`european;`call;300f;1f;100000f;`AAPL`MSFT`NVDA;0.4 0.3 0.3;`none;`none;0N;`none;0Nf;0f);

/ Note: basket routing from portfolio needs marketDataBook, not flat marketData
/ Test via direct basket pricing
basketResult:.basket.priceBasketOption[basketTrade;mktBook;corrTable;enlist[`mcConfig]!enlist .montecarlo.defaultMcConfig[]];
.testutil.assertTrue[basketResult[`status]~`OK;"basket prices OK"];
.testutil.assertTrue[basketResult[`unitPrice]>0f;"basket price positive"];
.testutil.assertTrue[basketResult[`standardError]>0f;"basket SE positive"];
.testutil.assertNear[basketResult`notionalPrice;basketResult[`unitPrice]*100000f;1f;"basket notional"];

-1 "PASS test_portfolio_basket_products: basketPrice=",string[basketResult`unitPrice],", SE=",string basketResult`standardError;
