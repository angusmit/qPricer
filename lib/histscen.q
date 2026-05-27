/ histscen.q - historical scenario shock table management (v0.28)

.histscen.validateHistoricalShockTable:{[historicalShockTable]
    if[0=count historicalShockTable; '"Historical shock table is empty"];
    firstRow:historicalShockTable 0;
    tableCols:key firstRow;
    requiredCols:`scenarioDate`eventName`underlying`spotReturn`volatilityShift`rateShift`dividendShift;
    missingCols:requiredCols where not requiredCols in tableCols;
    if[0<count missingCols; '"Missing columns: ",", " sv string missingCols];
    spotReturns:historicalShockTable`spotReturn;
    if[any spotReturns<= -1f; '"spotReturn must be > -1"];
 };

.histscen.scenarioKeys:{[historicalShockTable]
    allDates:historicalShockTable`scenarioDate;
    allEvents:historicalShockTable`eventName;
    pairKeys:();
    seen:();
    pIdx:0;
    while[pIdx<count allDates;
        pairKey:(allDates pIdx;allEvents pIdx);
        isNew:$[0=count seen;1b;not pairKey in seen];
        if[isNew;
            seen:seen,enlist pairKey;
            pairKeys:pairKeys,enlist `scenarioDate`eventName!(allDates pIdx;allEvents pIdx)];
        pIdx+:1];
    pairKeys
 };

.histscen.shocksForScenario:{[historicalShockTable;scenarioDate;eventName]
    dateMask:(historicalShockTable`scenarioDate)=scenarioDate;
    eventMask:(historicalShockTable`eventName)=eventName;
    bothMask:dateMask and eventMask;
    historicalShockTable where bothMask
 };

.histscen.applyHistoricalShock:{[marketDataBook;shockRows;configDict]
    missingMode:$[0<count configDict;
        $[`missingMode in key configDict;configDict`missingMode;`zeroShock];
        `zeroShock];
    / Copy market data book tables
    origSpotTable:marketDataBook`spotTable;
    origVolTable:marketDataBook`volatilityTable;
    origRateTable:marketDataBook`rateTable;
    origDivTable:marketDataBook`dividendTable;
    newSpots:origSpotTable`spot;
    newVols:origVolTable`volatility;
    newRates:origRateTable`riskFreeRate;
    newDivs:origDivTable`dividendYield;
    / Apply rate shock once (rates are global, not per-underlying)
    if[0<count shockRows;
        firstShock:shockRows 0;
        rateShiftVal:firstShock`rateShift;
        rateIdx:0;
        while[rateIdx<count newRates;
            newRates[rateIdx]:newRates[rateIdx]+rateShiftVal;
            rateIdx+:1]];
    / Apply per-underlying shocks
    shockIdx:0;
    while[shockIdx<count shockRows;
        shockRow:shockRows shockIdx;
        symbolName:shockRow`underlying;
        / Apply spot shock
        spotIdx:(origSpotTable`underlying)?symbolName;
        if[spotIdx<count origSpotTable;
            newSpots[spotIdx]:newSpots[spotIdx]*1f+shockRow`spotReturn];
        / Apply vol shock
        volIdx:(origVolTable`underlying)?symbolName;
        if[volIdx<count origVolTable;
            newVols[volIdx]:0.001|newVols[volIdx]+shockRow`volatilityShift];
        / Apply div shock
        divIdx:(origDivTable`underlying)?symbolName;
        if[divIdx<count origDivTable;
            newDivs[divIdx]:0f|newDivs[divIdx]+shockRow`dividendShift];
        shockIdx+:1];
    shockedSpotTable:([] underlying:origSpotTable`underlying; spot:newSpots);
    shockedVolTable:([] underlying:origVolTable`underlying; volatility:newVols);
    shockedRateTable:([] expiry:origRateTable`expiry; riskFreeRate:newRates);
    shockedDivTable:([] underlying:origDivTable`underlying; dividendYield:newDivs);
    shockedBook:.marketbook.createMarketDataBook[shockedSpotTable;shockedVolTable;shockedRateTable;shockedDivTable];
    shockedBook
 };

.histscen.applyAllHistoricalShocks:{[marketDataBook;historicalShockTable;configDict]
    scenKeys:.histscen.scenarioKeys historicalShockTable;
    resultRows:();
    sIdx:0;
    while[sIdx<count scenKeys;
        scenKey:scenKeys sIdx;
        shockRows:.histscen.shocksForScenario[historicalShockTable;scenKey`scenarioDate;scenKey`eventName];
        shockedBook:.histscen.applyHistoricalShock[marketDataBook;shockRows;configDict];
        resultRows:resultRows,enlist `scenarioDate`eventName`shockedMarketDataBook`status`errorMessage!(
            scenKey`scenarioDate;scenKey`eventName;shockedBook;`OK;"");
        sIdx+:1];
    resultRows
 };

.histscen.syntheticShockTable:{[symbolList]
    resultRows:();
    symIdx:0;
    while[symIdx<count symbolList;
        sym:symbolList symIdx;
        resultRows:resultRows,enlist `scenarioDate`eventName`underlying`spotReturn`volatilityShift`rateShift`dividendShift!(
            2020.03.16;`covidCrash;sym;-0.12;0.15;-0.005;0f);
        resultRows:resultRows,enlist `scenarioDate`eventName`underlying`spotReturn`volatilityShift`rateShift`dividendShift!(
            2022.06.13;`inflationShock;sym;-0.05;0.06;0.01;0f);
        resultRows:resultRows,enlist `scenarioDate`eventName`underlying`spotReturn`volatilityShift`rateShift`dividendShift!(
            2023.10.27;`rateHike;sym;-0.03;0.03;0.005;0f);
        symIdx+:1];
    resultRows
 };
