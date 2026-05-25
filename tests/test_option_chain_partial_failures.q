/ test_option_chain_partial_failures.q - bad rows do not crash the chain
\l lib/init.q

validCallPrice:.validation.blackScholesClosedForm[`call;100f;100f;1f;0.05;0f;0.2];
validPutPrice:.validation.blackScholesClosedForm[`put;100f;100f;1f;0.05;0f;0.2];

optionChain:([]
    underlying:`AAPL`AAPL`AAPL`AAPL;
    expiry:1 1 1 1f;
    strike:100 100 100 100f;
    optionType:`call`put`call`put;
    marketPrice:(validCallPrice;validPutPrice;0.0001;200f));

chainWithIVs:.iv.calculateOptionChainImpliedVols[optionChain;100f;0.05;0f];

/ 1. Row count preserved
if[not 4=count chainWithIVs; '"FAIL: expected 4 rows"];

/ 2. Count OK and ERROR rows
ivStatusList:chainWithIVs[;`ivStatus];
okRowCount:sum ivStatusList=`OK;
errorRowCount:sum ivStatusList=`ERROR;

if[not okRowCount=2; '"FAIL: expected 2 OK rows, got ",string okRowCount];
if[not errorRowCount=2; '"FAIL: expected 2 ERROR rows, got ",string errorRowCount];

/ 3. OK rows have positive implied vol
okIVs:chainWithIVs[;`impliedVolatility] where ivStatusList=`OK;
if[any okIVs<=0f; '"FAIL: OK rows should have positive implied vol"];

/ 4. ERROR rows have null implied vol
errorIVs:chainWithIVs[;`impliedVolatility] where ivStatusList=`ERROR;
if[not all null errorIVs; '"FAIL: ERROR rows should have null implied vol"];

/ 5. ERROR rows have non-empty status message
errorMsgs:chainWithIVs[;`ivStatusMessage] where ivStatusList=`ERROR;
if[any {0=count x} each errorMsgs; '"FAIL: ERROR rows should have non-empty status message"];

-1 "PASS test_option_chain_partial_failures: okRows=",string[okRowCount],", errorRows=",string errorRowCount;
