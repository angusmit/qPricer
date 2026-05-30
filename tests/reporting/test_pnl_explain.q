/ test_pnl_explain.q - PnL explain via Taylor expansion
\l core/init.q

/ European call and put - Greeks supported
tradeTable:([]
    tradeId:1 2;
    underlying:`AAPL`AAPL;
    productType:`equityOption`equityOption;
    exerciseStyle:`european`european;
    optionType:`call`put;
    strike:100 100f;
    expiry:1 1f;
    notional:1 1f;
    barrierType:`none`none;
    barrierLevel:0N 0N;
    rebate:0 0f);

/ t0 market data book
spotTable0:([] underlying:enlist `AAPL; spot:enlist 100f);
volTable0:([] underlying:enlist `AAPL; volatility:enlist 0.20);
rateTable0:([] expiry:enlist 1f; riskFreeRate:enlist 0.05);
divTable0:([] underlying:enlist `AAPL; dividendYield:enlist 0f);
book0:.marketbook.createMarketDataBook[spotTable0;volTable0;rateTable0;divTable0];

/ t1 market data book: spot up 1%, vol up 0.5pp, rate unchanged
spotTable1:([] underlying:enlist `AAPL; spot:enlist 101f);
volTable1:([] underlying:enlist `AAPL; volatility:enlist 0.205);
rateTable1:([] expiry:enlist 1f; riskFreeRate:enlist 0.05);
divTable1:([] underlying:enlist `AAPL; dividendYield:enlist 0f);
book1:.marketbook.createMarketDataBook[spotTable1;volTable1;rateTable1;divTable1];

bsModel:.model.createBlackScholesModel[];
fdmConfig:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

pnlConfig:`model`fdmConfig`timeStepYears`bookName!(bsModel;fdmConfig;1%252;"testBook");

explainResult:.pnl.explainPortfolio[tradeTable;book0;book1;pnlConfig];

/ 1. Two rows
if[not 2=count explainResult; '"FAIL: expected 2 rows"];

/ 2. Both status OK
statusList:explainResult`status;
if[not all statusList=`OK; '"FAIL: not all status OK"];

/ 3. actualPnL = pv1 - pv0
pv0List:explainResult`pv0;
pv1List:explainResult`pv1;
actualList:explainResult`actualPnL;
computedActual:pv1List-pv0List;
if[any (abs actualList-computedActual)>1e-10; '"FAIL: actualPnL != pv1 - pv0"];

/ 4. explainedPnL = sum of components
deltaList:explainResult`deltaPnL;
gammaList:explainResult`gammaPnL;
vegaList:explainResult`vegaPnL;
rhoList:explainResult`rhoPnL;
thetaList:explainResult`thetaPnL;
explainedList:explainResult`explainedPnL;
computedExplained:deltaList+gammaList+vegaList+rhoList+thetaList;
if[any (abs explainedList-computedExplained)>1e-10; '"FAIL: explainedPnL != sum of components"];

/ 5. unexplainedPnL = actual - explained
unexplainedList:explainResult`unexplainedPnL;
computedUnexplained:actualList-explainedList;
if[any (abs unexplainedList-computedUnexplained)>1e-10; '"FAIL: unexplainedPnL != actual - explained"];

/ 6. For small moves, unexplained should be small relative to actual
callUnexplained:abs unexplainedList 0;
callActual:abs actualList 0;
if[callActual>0.01;
    if[callUnexplained>0.5*callActual; '"FAIL: unexplained too large vs actual for call"]];

/ 7. Aggregate
aggResult:.pnl.aggregateExplain explainResult;
if[not aggResult[`totalActualPnL]=(sum actualList); '"FAIL: aggregate totalActualPnL mismatch"];

/ 8. Call delta PnL should be positive (spot went up, call delta > 0)
if[not deltaList[0]>0f; '"FAIL: call deltaPnL should be positive for spot up"];

/ 9. Put delta PnL should be negative (spot went up, put delta < 0)
if[not deltaList[1]<0f; '"FAIL: put deltaPnL should be negative for spot up"];

-1 "PASS test_pnl_explain: rows=2, callActualPnL=",string[actualList 0],", callUnexplained=",string[unexplainedList 0],", putActualPnL=",string actualList 1;
