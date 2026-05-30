/ test_surface_market_data_pricing.q - price using vol surface market data
\l core/init.q

/ Build surface where ATM 1Y vol = 0.2
spot:100f;
riskFreeRate:0.05;
dividendYield:0f;

/ Simple chain: all vol = 0.2
strikes:90 100 110f;
expiries:1 1 1f;
optTypes:`call`call`call;
chainPrices:{[ot;stk;expiryVal] .validation.blackScholesClosedForm[ot;spot;stk;expiryVal;riskFreeRate;dividendYield;0.2]}'[optTypes;strikes;expiries];

optionChain:([]
    underlying:`AAPL`AAPL`AAPL;
    expiry:expiries;
    strike:strikes;
    optionType:optTypes;
    marketPrice:chainPrices);

chainWithIVs:.iv.calculateOptionChainImpliedVols[optionChain;spot;riskFreeRate;dividendYield];
volSurface:.surface.buildVolSurface[chainWithIVs];
surfaceMarketData:.surface.createSurfaceMarketData[`AAPL;spot;riskFreeRate;dividendYield;volSurface];

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1f);
model:.model.createBlackScholesModel[];
config:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

surfaceResult:.engine.priceOption[trade;surfaceMarketData;model;config];
surfacePrice:surfaceResult`unitPrice;

/ Compare with flat market data
flatMkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;spot;riskFreeRate;dividendYield;0.2);
flatResult:.engine.priceOption[trade;flatMkt;model;config];
flatPrice:flatResult`unitPrice;

priceDiff:abs surfacePrice-flatPrice;
if[not surfacePrice>0f; '"FAIL: surface price not positive"];
if[priceDiff>0.01; '"FAIL: surface vs flat price diff too large: ",string priceDiff];

-1 "PASS test_surface_market_data_pricing: surfacePrice=",string[surfacePrice],", flatPrice=",string[flatPrice],", diff=",string priceDiff;
