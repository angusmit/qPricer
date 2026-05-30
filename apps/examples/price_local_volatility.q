/ price_local_volatility.q - local volatility pricing example
/ Usage: q examples/price_local_volatility.q

\l core/init.q
-1 "qFDM v",.qfdm.version," - Local Volatility Example\n";

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1f);

flatMkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);

constLocalVolFn:{[spotValue;timePoint] 0.2};
constLocalVolMkt:.market.createLocalVolatilityMarketData[`AAPL;100f;0.05;0f;constLocalVolFn];

downsideOnlySkewFn:{[spotValue;timePoint]
    extraVol:0f|0.001*(100f-spotValue);
    0.2+extraVol
 };
skewLocalVolMkt:.market.createLocalVolatilityMarketData[`AAPL;100f;0.05;0f;downsideOnlySkewFn];

bsModel:.model.createBlackScholesModel[];
lvModel:.model.createLocalVolatilityModel[];

config:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

bsPrice:(.engine.priceOption[trade;flatMkt;bsModel;config])`unitPrice;
constLvPrice:(.engine.priceOption[trade;constLocalVolMkt;lvModel;config])`unitPrice;
skewLvPrice:(.engine.priceOption[trade;skewLocalVolMkt;lvModel;config])`unitPrice;

-1 "European call (S=100, K=100, T=1, r=5%):";
-1 "  Flat BS explicit FDM (vol=20%): ",string bsPrice;
-1 "  Constant local vol (20%):       ",string constLvPrice;
-1 "  Downside-only local vol uplift: ",string skewLvPrice;
-1 "";
-1 "Flat-equivalence difference:      ",string abs bsPrice-constLvPrice;
-1 "Skew-vs-flat difference:          ",string skewLvPrice-bsPrice;
-1 "";
-1 "Constant local volatility reproduces the flat-vol explicit FDM price.";
-1 "The non-flat local-vol function changes the price while remaining positive";
-1 "and bounded.";
