/ test_vol_surface_lookup.q - build surface and test nearest-neighbour lookup
\l core/init.q

/ Build chain with varying vols
spot:100f;
riskFreeRate:0.05;
dividendYield:0f;

/ Create chain with different true vols per row
trueVols:0.22 0.20 0.21 0.23;
strikes:90 100 100 110f;
expiries:0.5 0.5 1 1f;
optTypes:`call`call`put`put;
chainPrices:{[ot;stk;expiryVal;vol] .validation.blackScholesClosedForm[ot;spot;stk;expiryVal;riskFreeRate;dividendYield;vol]}'[optTypes;strikes;expiries;trueVols];

optionChain:([]
    underlying:`AAPL`AAPL`AAPL`AAPL;
    expiry:expiries;
    strike:strikes;
    optionType:optTypes;
    marketPrice:chainPrices);

chainWithIVs:.iv.calculateOptionChainImpliedVols[optionChain;spot;riskFreeRate;dividendYield];
volSurface:.surface.buildVolSurface[chainWithIVs];

/ 1. Surface not empty
if[0=count volSurface; '"FAIL: vol surface is empty"];

/ 2. Exact lookup: strike=100, expiry=0.5
exactVol:.surface.getSurfaceVolatility[volSurface;100f;0.5];
if[(abs exactVol-0.20)>0.001; '"FAIL: exact lookup vol unexpected: ",string exactVol];

/ 3. Nearest lookup: strike=102, expiry=0.5 should match strike=100
nearVol:.surface.getSurfaceVolatility[volSurface;102f;0.5];
if[(abs nearVol-0.20)>0.001; '"FAIL: nearest lookup vol unexpected: ",string nearVol];

/ 4. Nearest lookup: strike=100, expiry=0.8 should match expiry=0.5 or 1.0
midExpiryVol:.surface.getSurfaceVolatility[volSurface;100f;0.8];
if[not midExpiryVol>0f; '"FAIL: mid-expiry lookup returned non-positive vol"];

-1 "PASS test_vol_surface_lookup: exactVol=",string[exactVol],", nearVol=",string nearVol;
