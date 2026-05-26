/ test_batch.q - daily batch orchestration
\l lib/init.q

tradeTable:([]
    tradeId:1 2;
    underlying:`AAPL`AAPL;
    productType:`equityOption`equityOption;
    exerciseStyle:`european`european;
    optionType:`call`put;
    strike:100 100f;
    expiry:1 1f;
    notional:1000000 500000f;
    barrierType:`none`none;
    barrierLevel:0N 0N;
    rebate:0 0f);

spotT0:([] underlying:enlist `AAPL; spot:enlist 100f);
volT0:([] underlying:enlist `AAPL; volatility:enlist 0.20);
rateT0:([] expiry:enlist 1f; riskFreeRate:enlist 0.05);
divT0:([] underlying:enlist `AAPL; dividendYield:enlist 0f);
book0:.marketbook.createMarketDataBook[spotT0;volT0;rateT0;divT0];

spotT1:([] underlying:enlist `AAPL; spot:enlist 101f);
volT1:([] underlying:enlist `AAPL; volatility:enlist 0.205);
rateT1:([] expiry:enlist 1f; riskFreeRate:enlist 0.05);
divT1:([] underlying:enlist `AAPL; dividendYield:enlist 0f);
book1:.marketbook.createMarketDataBook[spotT1;volT1;rateT1;divT1];

bsModel:.model.createBlackScholesModel[];
fdmConfig:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

configDict:`model`fdmConfig`timeStepYears`bookName`valuationDate`runLabel!(
    bsModel;fdmConfig;1%252;"testBook";2025.01.02;"dailyRun");

runResult:.batch.runDailyPricing[tradeTable;book1;book0;configDict];

/ 1. All expected keys present
expectedKeys:`pricingResult`greekResult`scenarioResult`pnlExplainResult`portfolioSummary`riskSummary`scenarioSummary`errorSummary`auditRecord;
resultKeys:key runResult;
missingKeys:expectedKeys where not expectedKeys in resultKeys;
if[0<count missingKeys; '"FAIL: missing keys: ",", " sv string missingKeys];

/ 2. Pricing has 2 rows
if[not 2=count runResult`pricingResult; '"FAIL: expected 2 pricing rows"];

/ 3. Scenarios: 2 trades * 9 = 18
if[not 18=count runResult`scenarioResult; '"FAIL: expected 18 scenario rows"];

/ 4. PnL explain has 2 rows
if[not 2=count runResult`pnlExplainResult; '"FAIL: expected 2 PnL rows"];

/ 5. Audit record has correct trade count
auditRec:runResult`auditRecord;
if[not auditRec[`tradeCount]=2; '"FAIL: audit tradeCount should be 2"];
if[not auditRec[`pricedOkCount]=2; '"FAIL: audit pricedOkCount should be 2"];

/ 6. Write reports
.batch.writeDailyReports[runResult;"/home/claude";"2025.01.02"];

-1 "PASS test_batch: keys=",string[count resultKeys],", pricing=2, scenarios=18, pnl=2";
