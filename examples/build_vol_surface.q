/ build_vol_surface.q - build and query a vol surface
/ Usage: q examples/build_vol_surface.q

\l lib/init.q
-1 "qFDM v",.qfdm.version," - Vol Surface Example\n";

spot:100f;
riskFreeRate:0.05;
dividendYield:0f;

/ Build option chain with varying true vols (simple skew)
trueVols:0.25 0.22 0.20 0.19 0.22 0.21 0.20 0.21;
strikes:80 90 100 110 80 90 100 110f;
expiries:0.5 0.5 0.5 0.5 1 1 1 1f;
optTypes:`call`call`call`call`put`put`put`put;
chainPrices:{[ot;stk;exp;vol] .validation.blackScholesClosedForm[ot;spot;stk;exp;riskFreeRate;dividendYield;vol]}'[optTypes;strikes;expiries;trueVols];

optionChain:([]
    underlying:8#`AAPL;
    expiry:expiries;
    strike:strikes;
    optionType:optTypes;
    marketPrice:chainPrices);

-1 "Option chain:";
show optionChain;
-1 "";

chainWithIVs:.iv.calculateOptionChainImpliedVols[optionChain;spot;riskFreeRate;dividendYield];

-1 "Implied volatilities:";
-1 "  "," " sv {string[x`strike],"@",string[x`expiry],"=",string x`impliedVolatility} each chainWithIVs;
-1 "";

volSurface:.surface.buildVolSurface[chainWithIVs];
-1 "Vol surface:";
show volSurface;
-1 "";

/ Query surface
surfaceVol100:.surface.getSurfaceVolatility[volSurface;100f;1f];
surfaceVol95:.surface.getSurfaceVolatility[volSurface;95f;0.5];
-1 "Surface vol at K=100, T=1: ",string surfaceVol100;
-1 "Surface vol at K=95, T=0.5 (nearest): ",string surfaceVol95;
-1 "";

/ Price using surface market data
surfaceMkt:.surface.createSurfaceMarketData[`AAPL;spot;riskFreeRate;dividendYield;volSurface];
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1f);
model:.model.createBlackScholesModel[];
config:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

surfacePrice:(.engine.priceOption[trade;surfaceMkt;model;config])`unitPrice;
-1 "European call priced from vol surface: ",string surfacePrice;
