/ marketbook.q - multi-symbol market data book
/ Routes each trade to the correct spot/vol/rate/dividend by underlying.

/ --- Public ---

.marketbook.createMarketDataBook:{[spotTable;volatilityTable;rateTable;dividendTable]
    marketDataBook:`spotTable`volatilityTable`rateTable`dividendTable!(spotTable;volatilityTable;rateTable;dividendTable);
    .marketbook.validateMarketDataBook marketDataBook;
    marketDataBook
 };

.marketbook.validateMarketDataBook:{[marketDataBook]
    .marketbook.__validateSpotTable marketDataBook`spotTable;
    .marketbook.__validateVolatilityTable marketDataBook`volatilityTable;
    .marketbook.__validateRateTable marketDataBook`rateTable;
    .marketbook.__validateDividendTable marketDataBook`dividendTable;
 };

.marketbook.getSpot:{[marketDataBook;underlying]
    spotRow:.marketbook.__getSingleRowForUnderlying[marketDataBook`spotTable;underlying;"spotTable"];
    spotRow`spot
 };

.marketbook.getVolatility:{[marketDataBook;underlying;strike;expiry]
    volRow:.marketbook.__getSingleRowForUnderlying[marketDataBook`volatilityTable;underlying;"volatilityTable"];
    volRow`volatility
 };

.marketbook.getRiskFreeRate:{[marketDataBook;expiry]
    rateTable:marketDataBook`rateTable;
    if[0=count rateTable; '"rateTable is empty"];
    expiryDistances:abs (rateTable`expiry)-expiry;
    nearestDist:min expiryDistances;
    matchRows:rateTable where expiryDistances=nearestDist;
    nearestRow:matchRows 0;
    nearestRow`riskFreeRate
 };

.marketbook.getDividendYield:{[marketDataBook;underlying;expiry]
    divRow:.marketbook.__getSingleRowForUnderlying[marketDataBook`dividendTable;underlying;"dividendTable"];
    divRow`dividendYield
 };

.marketbook.getMarketDataForTrade:{[marketDataBook;trade]
    tradeUnderlying:trade`underlying;
    tradeStrike:trade`strike;
    tradeExpiry:trade`expiry;
    spotVal:.marketbook.getSpot[marketDataBook;tradeUnderlying];
    volVal:.marketbook.getVolatility[marketDataBook;tradeUnderlying;tradeStrike;tradeExpiry];
    rateVal:.marketbook.getRiskFreeRate[marketDataBook;tradeExpiry];
    divVal:.marketbook.getDividendYield[marketDataBook;tradeUnderlying;tradeExpiry];
    `underlying`spot`riskFreeRate`dividendYield`volatility!(tradeUnderlying;spotVal;rateVal;divVal;volVal)
 };

/ --- Internal ---

.marketbook.__requireTableColumns:{[tbl;requiredCols;tableName]
    tableCols:cols tbl;
    missingCols:requiredCols where not requiredCols in tableCols;
    if[0<count missingCols; '"Missing columns in ",tableName,": ",", " sv string missingCols];
 };

.marketbook.__getSingleRowForUnderlying:{[tbl;underlying;tableName]
    matchRows:tbl where tbl[`underlying]=underlying;
    if[0=count matchRows; '"Symbol not found in ",tableName,": ",string underlying];
    matchRows 0
 };

.marketbook.__validateSpotTable:{[spotTable]
    .marketbook.__requireTableColumns[spotTable;`underlying`spot;"spotTable"];
    if[any (spotTable`spot)<=0f; '"spotTable contains non-positive spot"];
 };

.marketbook.__validateVolatilityTable:{[volatilityTable]
    .marketbook.__requireTableColumns[volatilityTable;`underlying`volatility;"volatilityTable"];
    if[any (volatilityTable`volatility)<=0f; '"volatilityTable contains non-positive volatility"];
 };

.marketbook.__validateRateTable:{[rateTable]
    .marketbook.__requireTableColumns[rateTable;`expiry`riskFreeRate;"rateTable"];
    if[any (rateTable`expiry)<=0f; '"rateTable contains non-positive expiry"];
 };

.marketbook.__validateDividendTable:{[dividendTable]
    .marketbook.__requireTableColumns[dividendTable;`underlying`dividendYield;"dividendTable"];
    if[any (dividendTable`dividendYield)<0f; '"dividendTable contains negative dividendYield"];
 };
