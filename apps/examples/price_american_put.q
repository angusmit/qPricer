/ price_american_put.q — price American put and compare to European
/ Usage: q examples/price_american_put.q

\l core/init.q
-1 "qFDM v",.qfdm.version," — American Put Pricing\n";

/ European put for comparison
europeanTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`european;`put;100f;1f;1f);

/ American put
americanTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    2;`AAPL;`equityOption;`american;`put;100f;1f;1f);

marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(
    `AAPL;100f;0.05;0f;0.2);

model:.model.createBlackScholesModel[];
config:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

-1 "European put:";
europeanResult:.engine.priceOption[europeanTrade;marketData;model;config];
show europeanResult;

-1 "\nAmerican put:";
americanResult:.engine.priceOption[americanTrade;marketData;model;config];
show americanResult;

earlyExercisePremium:americanResult[`unitPrice] - europeanResult`unitPrice;
-1 "\nEarly exercise premium: ",string earlyExercisePremium;
-1 "American/European ratio: ",string americanResult[`unitPrice] % europeanResult`unitPrice;
