/ test_portfolio_greeks.q - portfolio Greeks with unsupported handling
\l lib/init.q

tradeTable:([]
    tradeId:1 2 3 4;
    underlying:`AAPL`AAPL`AAPL`AAPL;
    productType:`equityOption`equityOption`equityOption`equityOption;
    exerciseStyle:`european`european`american`european;
    optionType:`call`put`put`call;
    strike:100 100 100 100f;
    expiry:1 1 1 1f;
    notional:1000000 1000000 1000000 1000000f;
    barrierType:`none`none`none`upAndOut;
    barrierLevel:0N 0N 0N 130f;
    rebate:0 0 0 0f);

marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
model:.model.createBlackScholesModel[];
config:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

portfolioGreeksTable:.portfolio.calculatePortfolioGreeks[tradeTable;marketData;model;config];

/ 1. Correct row count
if[not 4=count portfolioGreeksTable; '"FAIL: expected 4 rows"];

/ 2. European call and put are OK
statusValues:portfolioGreeksTable`status;
if[not statusValues[0]~`OK; '"FAIL: European call status should be OK"];
if[not statusValues[1]~`OK; '"FAIL: European put status should be OK"];

/ 3. American put and barrier call are UNSUPPORTED
if[not statusValues[2]~`UNSUPPORTED; '"FAIL: American put status should be UNSUPPORTED"];
if[not statusValues[3]~`UNSUPPORTED; '"FAIL: barrier call status should be UNSUPPORTED"];

/ 4. European call delta > 0
euroCallDelta:(portfolioGreeksTable`delta) 0;
if[not euroCallDelta>0f; '"FAIL: European call delta should be positive"];

/ 5. European put delta < 0
euroPutDelta:(portfolioGreeksTable`delta) 1;
if[not euroPutDelta<0f; '"FAIL: European put delta should be negative"];

/ 6. Both gammas positive
euroCallGamma:(portfolioGreeksTable`gamma) 0;
euroPutGamma:(portfolioGreeksTable`gamma) 1;
if[not euroCallGamma>0f; '"FAIL: European call gamma should be positive"];
if[not euroPutGamma>0f; '"FAIL: European put gamma should be positive"];

supportedRowCount:sum statusValues=`OK;
unsupportedRowCount:sum statusValues=`UNSUPPORTED;
-1 "PASS test_portfolio_greeks: rows=4, supportedRows=",string[supportedRowCount],", unsupportedRows=",string unsupportedRowCount;
