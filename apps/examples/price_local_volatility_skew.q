/ price_local_volatility_skew.q - non-flat local volatility example
/ Usage: q examples/price_local_volatility_skew.q

\l core/init.q
-1 "qFDM v",.qfdm.version," - Local Volatility Skew Example\n";

callTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1f);
putTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    2;`AAPL;`equityOption;`european;`put;100f;1f;1f);

flatMkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);

downsideOnlySkewFn:{[spotValue;timePoint]
    extraVol:0f|0.001*(100f-spotValue);
    0.2+extraVol
 };
skewMkt:.market.createLocalVolatilityMarketData[`AAPL;100f;0.05;0f;downsideOnlySkewFn];

bsModel:.model.createBlackScholesModel[];
lvModel:.model.createLocalVolatilityModel[];
config:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

flatCallPrice:(.engine.priceOption[callTrade;flatMkt;bsModel;config])`unitPrice;
skewCallPrice:(.engine.priceOption[callTrade;skewMkt;lvModel;config])`unitPrice;
flatPutPrice:(.engine.priceOption[putTrade;flatMkt;bsModel;config])`unitPrice;
skewPutPrice:(.engine.priceOption[putTrade;skewMkt;lvModel;config])`unitPrice;

-1 "Downside-only local volatility uplift:";
-1 "  vol = 0.2 + max(0, 0.001 * (100 - S))";
-1 "";
-1 "Examples:";
-1 "  S=100: vol=20%";
-1 "  S=50:  vol=25%";
-1 "  S>100: vol=20%";
-1 "";
-1 "optionType  flatVolPrice  localVolPrice  priceDifference";
-1 "---------- ------------- -------------- ----------------";
-1 "call       ",("  " sv string each (flatCallPrice;skewCallPrice;skewCallPrice-flatCallPrice));
-1 "put        ",("  " sv string each (flatPutPrice;skewPutPrice;skewPutPrice-flatPutPrice));
-1 "";
-1 "The local-vol function is never below the flat 20% volatility, so both call";
-1 "and put prices increase versus the flat-vol benchmark. The same local-vol";
-1 "surface is used for both products, so the call-put relationship remains";
-1 "consistent under the shared model assumptions.";
