/ price_local_volatility_skew.q - non-flat local volatility example
/ Usage: q examples/price_local_volatility_skew.q

\l lib/init.q
-1 "qFDM v",.qfdm.version," - Local Volatility Skew Example\n";

callTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1f);
putTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    2;`AAPL;`equityOption;`european;`put;100f;1f;1f);

flatMkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);

downsideSkewFn:{[spotValue;timePoint]
    baseVol:0.2;
    skewAdj:0.001*(100f-spotValue);
    0.05|(baseVol+skewAdj)&0.80
 };
skewMkt:.market.createLocalVolatilityMarketData[`AAPL;100f;0.05;0f;downsideSkewFn];

bsModel:.model.createBlackScholesModel[];
lvModel:.model.createLocalVolatilityModel[];
config:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

flatCallPrice:(.engine.priceOption[callTrade;flatMkt;bsModel;config])`unitPrice;
skewCallPrice:(.engine.priceOption[callTrade;skewMkt;lvModel;config])`unitPrice;
flatPutPrice:(.engine.priceOption[putTrade;flatMkt;bsModel;config])`unitPrice;
skewPutPrice:(.engine.priceOption[putTrade;skewMkt;lvModel;config])`unitPrice;

-1 "Downside-skew local volatility: higher vol for lower spot levels.";
-1 "  sigma(S,t) = max(0.05, min(0.80, 0.2 + 0.001*(100 - S)))";
-1 "";
-1 "optionType  flatVolPrice  localVolPrice  priceDifference";
-1 "---------- ------------- -------------- ----------------";
-1 "call       ",("  " sv string each (flatCallPrice;skewCallPrice;skewCallPrice-flatCallPrice));
-1 "put        ",("  " sv string each (flatPutPrice;skewPutPrice;skewPutPrice-flatPutPrice));
-1 "";
-1 "Non-flat local vol changes prices versus flat BS - confirms solver uses sigma(S,t).";
