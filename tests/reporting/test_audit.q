/ test_audit.q - audit record correctness
\l core/init.q

tradeTable:([]
    tradeId:1 2 3;
    underlying:`AAPL`AAPL`TSLA;
    productType:`equityOption`equityOption`equityOption;
    exerciseStyle:`european`european`european;
    optionType:`call`put`call;
    strike:100 100 200f;
    expiry:1 1 1f;
    notional:1000000 500000 1000000f;
    barrierType:`none`none`none;
    barrierLevel:0N 0N 0N;
    rebate:0 0 0f);

/ Book only has AAPL - TSLA will fail
spotT0:([] underlying:enlist `AAPL; spot:enlist 100f);
volT0:([] underlying:enlist `AAPL; volatility:enlist 0.20);
rateT0:([] expiry:enlist 1f; riskFreeRate:enlist 0.05);
divT0:([] underlying:enlist `AAPL; dividendYield:enlist 0f);
book0:.marketbook.createMarketDataBook[spotT0;volT0;rateT0;divT0];
book1:.marketbook.createMarketDataBook[spotT0;volT0;rateT0;divT0];

bsModel:.model.createBlackScholesModel[];
fdmConfig:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

configDict:`model`fdmConfig`timeStepYears`bookName`valuationDate`runLabel!(
    bsModel;fdmConfig;1%252;"testBook";2025.01.02;"auditTest");

runResult:.batch.runDailyPricing[tradeTable;book1;book0;configDict];

/ 1. Audit counts
auditRec:runResult`auditRecord;
if[not auditRec[`tradeCount]=3; '"FAIL: tradeCount should be 3"];
if[not auditRec[`pricedOkCount]=2; '"FAIL: pricedOkCount should be 2 (AAPL trades)"];
if[not auditRec[`errorCount]=1; '"FAIL: errorCount should be 1 (TSLA)"];

/ 2. Error summary
errorSummary:runResult`errorSummary;
if[not errorSummary[`pricingErrorCount]=1; '"FAIL: pricing errors should be 1"];

/ 3. Scenario rows: 2 OK trades * 9 + 1 error trade * 1 = 19
if[not auditRec[`totalScenarioRows]=19; '"FAIL: totalScenarioRows should be 19"];

/ 4. PnL explain: 2 OK + 1 error = 3
if[not auditRec[`totalPnlExplainRows]=3; '"FAIL: totalPnlExplainRows should be 3"];

-1 "PASS test_audit: tradeCount=",string[auditRec`tradeCount],", okCount=",string[auditRec`pricedOkCount],", errorCount=",string auditRec`errorCount;
