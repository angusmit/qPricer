/ test_option_chain_implied_vols.q - implied vols for an option chain
\l core/init.q

trueVol:0.2;
spot:100f;
riskFreeRate:0.05;
dividendYield:0f;

/ Build chain with known BS prices
strikes:90 100 100 110f;
expiries:0.5 0.5 1 1f;
optTypes:`call`call`put`put;
chainPrices:{[ot;stk;expiryVal] .validation.blackScholesClosedForm[ot;spot;stk;expiryVal;riskFreeRate;dividendYield;trueVol]}'[optTypes;strikes;expiries];

optionChain:([]
    underlying:`AAPL`AAPL`AAPL`AAPL;
    expiry:expiries;
    strike:strikes;
    optionType:optTypes;
    marketPrice:chainPrices);

chainWithIVs:.iv.calculateOptionChainImpliedVols[optionChain;spot;riskFreeRate;dividendYield];

/ 1. Row count unchanged
if[not 4=count chainWithIVs; '"FAIL: expected 4 rows"];

/ 2. All IV status OK
ivStatusValues:chainWithIVs[;`ivStatus];
if[not all ivStatusValues=`OK; '"FAIL: not all ivStatus OK"];

/ 3. All implied vols close to 0.2
ivValues:chainWithIVs[;`impliedVolatility];
ivErrors:{abs x-trueVol} each ivValues;
if[any ivErrors>1e-4; '"FAIL: some implied vols too far from true vol"];

avgIV:avg ivValues;
-1 "PASS test_option_chain_implied_vols: rows=4, averageIV=",string avgIV;
