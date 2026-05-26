/ basket.q - basket option pricing via correlated Monte Carlo (v0.17)

/ --- Public ---

.basket.priceBasketOption:{[trade;marketDataBook;correlationTable;configDict]
    .basket.validateBasketTrade trade;
    basketSymbols:trade`basketSymbols;
    basketWeights:trade`basketWeights;
    numSymbols:count basketSymbols;
    mcConfig:$[`mcConfig in key configDict; configDict`mcConfig; .montecarlo.defaultMcConfig[]];
    .montecarlo.validateMcConfig mcConfig;
    / Extract market data vectors
    spotVector:numSymbols#0f;
    rateVector:numSymbols#0f;
    divVector:numSymbols#0f;
    volVector:numSymbols#0f;
    symIdx:0;
    while[symIdx<numSymbols;
          spotVector[symIdx]:.marketbook.getSpot[marketDataBook;basketSymbols symIdx];
          volVector[symIdx]:.marketbook.getVolatility[marketDataBook;basketSymbols symIdx;trade`strike;trade`expiry];
          rateVector[symIdx]:.marketbook.getRiskFreeRate[marketDataBook;trade`expiry];
          divVector[symIdx]:.marketbook.getDividendYield[marketDataBook;basketSymbols symIdx;trade`expiry];
          symIdx+:1];
    / Build correlation matrix
    .correlation.validateCorrelationTable[correlationTable;basketSymbols];
    correlationMatrix:.correlation.toMatrix[correlationTable;basketSymbols];
    / Simulate correlated paths
    pathData:.montecarlo.simulateCorrelatedGBMPaths[spotVector;rateVector;divVector;volVector;trade`expiry;correlationMatrix;mcConfig];
    / Terminal spots per symbol
    terminalSpots:{last each x} each pathData;
    / Basket value per path: sum(weight_i * S_i(T))
    weightedTerminals:basketWeights*terminalSpots;
    basketValues:sum weightedTerminals;
    / Payoff
    payoffVector:.basket.payoff[basketValues;trade`strike;trade`optionType];
    / Price
    avgRate:avg rateVector;
    priceResult:.montecarlo.priceFromPayoffs[payoffVector;avgRate;trade`expiry;mcConfig`confidenceLevel];
    `tradeId`underlying`optionType`productType`unitPrice`notionalPrice`method`modelName`standardError`lowerConfidence`upperConfidence`status`statusMessage!(
        trade`tradeId;trade[`basketSymbols]0;trade`optionType;`basketOption;
        priceResult`price;priceResult[`price]*trade`notional;
        `monteCarlo;`basketOption;
        priceResult`standardError;priceResult`lowerConfidence;priceResult`upperConfidence;
        `OK;"")
 };

.basket.basketValue:{[pathData;basketWeights]
    terminalSpots:{last each x} each pathData;
    sum basketWeights*terminalSpots
 };

.basket.payoff:{[basketValueVector;strike;optionType]
    if[optionType~`call; :0f|basketValueVector-strike];
    if[optionType~`put; :0f|strike-basketValueVector];
    '"Unsupported optionType: ",string optionType
 };

/ --- Greeks via bump-and-reprice ---

.basket.bumpGreek:{[trade;marketDataBook;correlationTable;configDict;greekName;targetSymbol]
    baseResult:.basket.priceBasketOption[trade;marketDataBook;correlationTable;configDict];
    basePrice:baseResult`unitPrice;
    bumpedBook:marketDataBook;
    bumpSize:0f;
    if[greekName~`delta;
       spotVal:.marketbook.getSpot[marketDataBook;targetSymbol];
       bumpSize:spotVal*0.01;
       origSpotTable:marketDataBook`spotTable;
       newSpots:origSpotTable`spot;
       spotIdx:(origSpotTable`underlying)?targetSymbol;
       newSpots[spotIdx]:spotVal+bumpSize;
       bumpedSpotTable:([] underlying:origSpotTable`underlying; spot:newSpots);
       bumpedBook:@[marketDataBook;`spotTable;:;bumpedSpotTable]];
    if[greekName~`vega;
       bumpSize:0.01;
       origVolTable:marketDataBook`volatilityTable;
       newVols:origVolTable`volatility;
       volIdx:(origVolTable`underlying)?targetSymbol;
       newVols[volIdx]:newVols[volIdx]+bumpSize;
       bumpedVolTable:([] underlying:origVolTable`underlying; volatility:newVols);
       bumpedBook:@[marketDataBook;`volatilityTable;:;bumpedVolTable]];
    if[bumpSize=0f; '"Unsupported greekName: ",string greekName];
    bumpedResult:.basket.priceBasketOption[trade;bumpedBook;correlationTable;configDict];
    (bumpedResult[`unitPrice]-basePrice)%bumpSize
 };

/ --- Validation ---

.basket.validateBasketTrade:{[trade]
    if[not trade[`productType]~`basketOption; '"productType must be `basketOption"];
    if[not trade[`optionType] in `call`put; '"optionType must be `call or `put"];
    if[not trade[`exerciseStyle]~`european; '"Basket options only support European exercise"];
    if[not trade[`strike]>0f; '"strike must be positive"];
    if[not trade[`expiry]>0f; '"expiry must be positive"];
    basketSymbols:trade`basketSymbols;
    basketWeights:trade`basketWeights;
    if[0=count basketSymbols; '"basketSymbols must be non-empty"];
    if[not (count basketSymbols)=count basketWeights; '"basketSymbols and basketWeights length mismatch"];
    if[(sum basketWeights)<=0f; '"sum of basketWeights must be positive"];
 };
