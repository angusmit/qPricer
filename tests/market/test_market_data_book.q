/ test_market_data_book.q - market data book creation and lookup
\l lib/init.q

spotTable:([] underlying:`AAPL`MSFT`NVDA; spot:100 250 800f);
volatilityTable:([] underlying:`AAPL`MSFT`NVDA; volatility:0.20 0.25 0.35);
rateTable:([] expiry:0.25 0.5 1 2f; riskFreeRate:0.045 0.0475 0.05 0.0525);
dividendTable:([] underlying:`AAPL`MSFT`NVDA; dividendYield:0 0.01 0f);

marketDataBook:.marketbook.createMarketDataBook[spotTable;volatilityTable;rateTable;dividendTable];

/ 1. Spot lookups
aaplSpot:.marketbook.getSpot[marketDataBook;`AAPL];
msftSpot:.marketbook.getSpot[marketDataBook;`MSFT];
if[not aaplSpot=100f; '"FAIL: AAPL spot expected 100"];
if[not msftSpot=250f; '"FAIL: MSFT spot expected 250"];

/ 2. Volatility lookup
nvdaVol:.marketbook.getVolatility[marketDataBook;`NVDA;800f;1f];
if[not nvdaVol=0.35; '"FAIL: NVDA vol expected 0.35"];

/ 3. Rate lookup (nearest expiry)
rate1Y:.marketbook.getRiskFreeRate[marketDataBook;1f];
if[not rate1Y=0.05; '"FAIL: 1Y rate expected 0.05"];

/ 4. Dividend lookup
msftDiv:.marketbook.getDividendYield[marketDataBook;`MSFT;1f];
if[not msftDiv=0.01; '"FAIL: MSFT div expected 0.01"];

/ 5. getMarketDataForTrade
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1f;`none;0Nf;0f);
marketDataForTrade:.marketbook.getMarketDataForTrade[marketDataBook;trade];
if[not marketDataForTrade[`underlying]~`AAPL; '"FAIL: underlying mismatch"];
if[not marketDataForTrade[`spot]=100f; '"FAIL: spot mismatch"];
if[not marketDataForTrade[`volatility]=0.20; '"FAIL: volatility mismatch"];
if[not marketDataForTrade[`riskFreeRate]=0.05; '"FAIL: rate mismatch"];
if[not marketDataForTrade[`dividendYield]=0f; '"FAIL: dividend mismatch"];

-1 "PASS test_market_data_book: AAPL spot=",string[aaplSpot],", MSFT spot=",string[msftSpot],", NVDA vol=",string nvdaVol;
