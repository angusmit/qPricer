/ test_basket_greeks.q - basket Greeks via bump-and-reprice
\l lib/init.q

spotTable:([] underlying:`AAPL`MSFT; spot:100 200f);
volTable:([] underlying:`AAPL`MSFT; volatility:0.2 0.25);
rateTable:([] expiry:enlist 1f; riskFreeRate:enlist 0.05);
divTable:([] underlying:`AAPL`MSFT; dividendYield:0 0.01f);
mktBook:.marketbook.createMarketDataBook[spotTable;volTable;rateTable;divTable];
corrTable:([] sym1:enlist `AAPL; sym2:enlist `MSFT; correlation:enlist 0.5);

callTrade:`tradeId`productType`basketSymbols`basketWeights`optionType`exerciseStyle`strike`expiry`notional!(
    1;`basketOption;`AAPL`MSFT;0.5 0.5;`call;`european;150f;1f;1f);

mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(
    50000;10;42;0b;0b;0.95);
configDict:enlist[`mcConfig]!enlist mcConfig;

/ 1. Call delta for AAPL should be positive
aaplDelta:.basket.bumpGreek[callTrade;mktBook;corrTable;configDict;`delta;`AAPL];
.testutil.assertTrue[aaplDelta>0f;"AAPL call delta positive"];

/ 2. Call delta for MSFT should be positive
msftDelta:.basket.bumpGreek[callTrade;mktBook;corrTable;configDict;`delta;`MSFT];
.testutil.assertTrue[msftDelta>0f;"MSFT call delta positive"];

/ 3. Vega should be positive
aaplVega:.basket.bumpGreek[callTrade;mktBook;corrTable;configDict;`vega;`AAPL];
.testutil.assertTrue[aaplVega>0f;"AAPL vega positive"];

/ 4. Same seed makes Greeks stable (recompute)
aaplDelta2:.basket.bumpGreek[callTrade;mktBook;corrTable;configDict;`delta;`AAPL];
.testutil.assertNear[aaplDelta;aaplDelta2;0.001;"delta stable with same seed"];

-1 "PASS test_basket_greeks: AAPL delta=",string[aaplDelta],", MSFT delta=",string msftDelta;
