/ test_regression.q - regression warnings
\l core/init.q

/ Build two pricing runs with different prices
currentPricing:([]
    tradeId:1 2 3;
    underlying:`AAPL`AAPL`MSFT;
    unitPrice:10.5 5.6 30.2f;
    status:`OK`OK`OK);

previousPricing:([]
    tradeId:1 2;
    underlying:`AAPL`AAPL;
    unitPrice:10.45 5.55f;
    status:`OK`OK);

/ 1. Compare with small tolerance - should flag trade 2
regressionResult:.regression.comparePricingRuns[currentPricing;previousPricing;0.04];
if[not 3=count regressionResult; '"FAIL: expected 3 regression rows"];

/ Trade 1: change = 0.05, below tolerance 0.04? No, 0.05 > 0.04 -> WARNING
trade1Status:(regressionResult`status) 0;
trade2Status:(regressionResult`status) 1;
trade3Status:(regressionResult`status) 2;

/ Trade 1 change = 0.05 > 0.04 -> WARNING
if[not trade1Status~`WARNING; '"FAIL: trade 1 should be WARNING (change 0.05 > 0.04)"];
/ Trade 2 change = 0.05 > 0.04 -> WARNING
if[not trade2Status~`WARNING; '"FAIL: trade 2 should be WARNING"];
/ Trade 3 no previous -> NO_PREVIOUS
if[not trade3Status~`NO_PREVIOUS; '"FAIL: trade 3 should be NO_PREVIOUS"];

/ 2. Compare with large tolerance - no warnings
regressionResult2:.regression.comparePricingRuns[currentPricing;previousPricing;1.0];
warningCount:sum (regressionResult2`status)=`WARNING;
if[not warningCount=0; '"FAIL: large tolerance should produce no warnings"];

/ 3. Flag large PnL
pnlResult:([]
    tradeId:1 2;
    actualPnL:0.5 0.01f;
    unexplainedPnL:0.1 0.001f;
    status:`OK`OK);

largePnLRows:.regression.flagLargePnL[pnlResult;0.1];
if[not 1=count largePnLRows; '"FAIL: should flag 1 large PnL trade"];

/ 4. Flag large unexplained
largeUnexplainedRows:.regression.flagLargeUnexplainedPnL[pnlResult;0.05];
if[not 1=count largeUnexplainedRows; '"FAIL: should flag 1 large unexplained trade"];

-1 "PASS test_regression: comparisons=",string[count regressionResult],", largePnL=",string[count largePnLRows],", largeUnexplained=",string count largeUnexplainedRows;
