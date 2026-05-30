/ analyze_american_put.q - American put pricing and early exercise analysis
/ Usage: q examples/analyze_american_put.q

\l core/init.q
-1 "qFDM v",.qfdm.version," - American Put Early Exercise Analysis\n";

americanTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`american;`put;100f;1f;1f);
europeanTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    2;`AAPL;`equityOption;`european;`put;100f;1f;1f);

marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(
    `AAPL;100f;0.05;0f;0.2);
model:.model.createBlackScholesModel[];
config:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

-1 "American put price:";
americanResult:.engine.priceOption[americanTrade;marketData;model;config];
show americanResult;

-1 "\nEuropean put price:";
europeanResult:.engine.priceOption[europeanTrade;marketData;model;config];
show europeanResult;

earlyExercisePremium:americanResult[`unitPrice]-europeanResult`unitPrice;
-1 "\nEarly exercise premium: ",string earlyExercisePremium;

-1 "\nEarly exercise boundary (sampled):";
boundaryTable:.american.extractEarlyExerciseBoundary[americanTrade;marketData;model;config];
exerciseRows:boundaryTable where boundaryTable`hasExerciseRegion;
/ Show every 200th row for readability plus the last row
sampleIndices:200*til 1+(-1+count exerciseRows) div 200;
sampleIndices:distinct sampleIndices,(-1)+count exerciseRows;
show exerciseRows sampleIndices;

-1 "\nFull analysis:";
analysisResult:.american.analyzeAmericanPut[americanTrade;marketData;model;config];
show analysisResult`earlyExercisePremium;
