/ price_local_volatility.q - local volatility pricing example
/ Usage: q examples/price_local_volatility.q

\l lib/init.q
-1 "qFDM v",.qfdm.version," - Local Volatility Example\n";

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1f);

/ Flat BS market data
flatMkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);

/ Local vol (constant) - should match flat BS exactly
constLocalVolFn:{[spotValue;timePoint] 0.2};
constLocalVolMkt:.market.createLocalVolatilityMarketData[`AAPL;100f;0.05;0f;constLocalVolFn];

/ Local vol with skew
skewLocalVolFn:{[spotValue;timePoint]
    baseVol:0.2;
    skewAdj:0.0005*100f-spotValue;
    candidateVol:baseVol+skewAdj;
    0.05|candidateVol
 };
skewLocalVolMkt:.market.createLocalVolatilityMarketData[`AAPL;100f;0.05;0f;skewLocalVolFn];

bsModel:.model.createBlackScholesModel[];
lvModel:.model.createLocalVolatilityModel[];

config:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

bsPrice:(.engine.priceOption[trade;flatMkt;bsModel;config])`unitPrice;
constLvPrice:(.engine.priceOption[trade;constLocalVolMkt;lvModel;config])`unitPrice;
skewLvPrice:(.engine.priceOption[trade;skewLocalVolMkt;lvModel;config])`unitPrice;

-1 "European call (S=100, K=100, T=1, r=5%):";
-1 "  Flat BS (vol=20%):         ",string bsPrice;
-1 "  Constant local vol (20%):  ",string constLvPrice;
-1 "  Skew local vol:            ",string skewLvPrice;
-1 "";
-1 "Flat equivalence diff:       ",string abs bsPrice-constLvPrice;
-1 "Skew vs flat diff:           ",string skewLvPrice-bsPrice;
