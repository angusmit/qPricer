/ price_barrier_options.q - barrier option pricing examples
/ Usage: q examples/price_barrier_options.q

\l lib/init.q
-1 "qFDM v",.qfdm.version," - Barrier Option Examples\n";

marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(
    `AAPL;100f;0.05;0f;0.2);
model:.model.createBlackScholesModel[];
config:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;300;4000;0f;300f;`linear;1b;1b);

/ --- Up-and-out call ---
vanillaCall:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1f);
barrierCall:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
    101;`AAPL;`equityOption;`european;`call;100f;1f;1f;`upAndOut;130f;0f);

vanillaCallResult:.engine.priceOption[vanillaCall;marketData;model;config];
barrierCallResult:.engine.priceOption[barrierCall;marketData;model;config];

-1 "Up-and-out call (barrier=130):";
-1 "  Vanilla call price:  ",string vanillaCallResult`unitPrice;
-1 "  Barrier call price:  ",string barrierCallResult`unitPrice;
-1 "  Barrier discount:    ",string vanillaCallResult[`unitPrice]-barrierCallResult`unitPrice;
-1 "";

/ --- Down-and-out put ---
vanillaPut:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    2;`AAPL;`equityOption;`european;`put;100f;1f;1f);
barrierPut:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
    102;`AAPL;`equityOption;`european;`put;100f;1f;1f;`downAndOut;70f;0f);

vanillaPutResult:.engine.priceOption[vanillaPut;marketData;model;config];
barrierPutResult:.engine.priceOption[barrierPut;marketData;model;config];

-1 "Down-and-out put (barrier=70):";
-1 "  Vanilla put price:   ",string vanillaPutResult`unitPrice;
-1 "  Barrier put price:   ",string barrierPutResult`unitPrice;
-1 "  Barrier discount:    ",string vanillaPutResult[`unitPrice]-barrierPutResult`unitPrice;
