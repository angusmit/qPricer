/ replay.q - historical scenario replay engine (v0.28)

/ Use Crank-Nicolson for replay — unconditionally stable under shocked vols
.replay.__stableConfig:{[]
    cfg:.config.defaultPricingConfig[];
    cfg[`method]:`crankNicolson;
    cfg
 };

.replay.basePortfolioValue:{[tradeTable;marketDataBook;configDict]
    bsModel:.model.createBlackScholesModel[];
    fdmConfig:.replay.__stableConfig[];
    resultRows:();
    tIdx:0;
    while[tIdx<count tradeTable;
        tradeRow:tradeTable tIdx;
        engineWrapper:{
            mkt:.marketbook.getMarketDataForTrade[x 1;x 0];
            .engine.priceOption[x 0;mkt;x 2;x 3]};
        priceResult:@[engineWrapper;(tradeRow;marketDataBook;bsModel;fdmConfig);{x}];
        pvVal:$[10h=type priceResult;0Nf;priceResult[`unitPrice]*tradeRow`notional];
        resultRows:resultRows,enlist `tradeId`basePV!(tradeRow`tradeId;pvVal);
        tIdx+:1];
    resultRows
 };

.replay.replayOneScenario:{[tradeTable;marketDataBook;historicalShockTable;scenarioDate;eventName;configDict]
    shockRows:.histscen.shocksForScenario[historicalShockTable;scenarioDate;eventName];
    shockWrapper:{.histscen.applyHistoricalShock[x 0;x 1;x 2]};
    shockedBookResult:@[shockWrapper;(marketDataBook;shockRows;configDict);{x}];
    if[10h=type shockedBookResult;
        / Shock application failed - return error rows for all trades
        resultRows:();
        errIdx:0;
        while[errIdx<count tradeTable;
            tradeRow:tradeTable errIdx;
            bookName:$[`bookName in key tradeRow;tradeRow`bookName;`default];
            resultRows:resultRows,enlist `scenarioDate`eventName`tradeId`bookName`underlying`productType`basePV`shockedPV`historicalPnl`status`errorMessage!(
                scenarioDate;eventName;tradeRow`tradeId;bookName;tradeRow`underlying;tradeRow`productType;0Nf;0Nf;0Nf;`ERROR;shockedBookResult);
            errIdx+:1];
        :resultRows];
    shockedBook:shockedBookResult;
    bsModel:.model.createBlackScholesModel[];
    fdmConfig:.replay.__stableConfig[];
    basePVs:.replay.basePortfolioValue[tradeTable;marketDataBook;configDict];
    resultRows:();
    tIdx:0;
    while[tIdx<count tradeTable;
        tradeRow:tradeTable tIdx;
        bookName:$[`bookName in key tradeRow;tradeRow`bookName;`default];
        / Shocked PV
        shockedWrapper:{
            mkt:.marketbook.getMarketDataForTrade[x 1;x 0];
            .engine.priceOption[x 0;mkt;x 2;x 3]};
        shockedResult:@[shockedWrapper;(tradeRow;shockedBook;bsModel;fdmConfig);{x}];
        basePVRow:basePVs tIdx;
        basePVVal:basePVRow`basePV;
        shockedPVVal:$[10h=type shockedResult;0Nf;shockedResult[`unitPrice]*tradeRow`notional];
        histPnl:$[(not null basePVVal) and not null shockedPVVal;shockedPVVal-basePVVal;0Nf];
        statusVal:$[null histPnl;`ERROR;`OK];
        errMsg:$[10h=type shockedResult;shockedResult;""];
        resultRows:resultRows,enlist `scenarioDate`eventName`tradeId`bookName`underlying`productType`basePV`shockedPV`historicalPnl`status`errorMessage!(
            scenarioDate;eventName;tradeRow`tradeId;bookName;tradeRow`underlying;tradeRow`productType;
            basePVVal;shockedPVVal;histPnl;statusVal;errMsg);
        tIdx+:1];
    resultRows
 };

.replay.replayHistoricalScenarios:{[tradeTable;marketDataBook;historicalShockTable;configDict]
    .histscen.validateHistoricalShockTable historicalShockTable;
    scenKeys:.histscen.scenarioKeys historicalShockTable;
    allResults:();
    sIdx:0;
    while[sIdx<count scenKeys;
        scenKey:scenKeys sIdx;
        scenResult:.replay.replayOneScenario[tradeTable;marketDataBook;historicalShockTable;scenKey`scenarioDate;scenKey`eventName;configDict];
        allResults:allResults,scenResult;
        sIdx+:1];
    allResults
 };

.replay.historicalPnlDistribution:{[replayResult]
    scenDates:replayResult`scenarioDate;
    scenEvents:replayResult`eventName;
    seen:();
    resultRows:();
    rIdx:0;
    while[rIdx<count replayResult;
        sd:scenDates rIdx;
        en:scenEvents rIdx;
        pairKey:(sd;en);
        isNew:$[0=count seen;1b;not pairKey in seen];
        if[isNew;
            seen:seen,enlist pairKey;
            dateMask:scenDates=sd;
            eventMask:scenEvents=en;
            bothMask:dateMask and eventMask;
            scenRows:replayResult where bothMask;
            okMask:(scenRows`status)=`OK;
            okRows:scenRows where okMask;
            totalPnl:$[0<count okRows;sum okRows`historicalPnl;0Nf];
            lossVal:$[not null totalPnl;neg totalPnl;0Nf];
            resultRows:resultRows,enlist `scenarioDate`eventName`pnl`loss`status`errorMessage!(
                sd;en;totalPnl;lossVal;`OK;"")];
        rIdx+:1];
    resultRows
 };

.replay.historicalVarReport:{[replayResult;confidenceLevels]
    pnlDist:.replay.historicalPnlDistribution replayResult;
    pnlVec:pnlDist`pnl;
    .var.varReport[pnlVec;confidenceLevels]
 };

.replay.worstHistoricalEvents:{[replayResult;topCountVal]
    pnlDist:.replay.historicalPnlDistribution replayResult;
    pnlCol:pnlDist`pnl;
    sortIdx:iasc pnlCol;
    takeCount:topCountVal&count pnlDist;
    resultRows:();
    rankIdx:0;
    while[rankIdx<takeCount;
        origIdx:sortIdx rankIdx;
        rowData:pnlDist origIdx;
        resultRows:resultRows,enlist `rank`scenarioDate`eventName`pnl`loss`status`errorMessage!(
            rankIdx+1;rowData`scenarioDate;rowData`eventName;rowData`pnl;rowData`loss;`OK;"");
        rankIdx+:1];
    resultRows
 };

.replay.historicalRiskContribution:{[replayResult;groupColumn;confidenceLevel]
    pnlTable:();
    rIdx:0;
    while[rIdx<count replayResult;
        rowData:replayResult rIdx;
        scenName:`$string[rowData`scenarioDate],"_",string rowData`eventName;
        groupVal:$[groupColumn in key rowData;rowData groupColumn;`unknown];
        pnlTable:pnlTable,enlist `scenarioName`pnl!(scenName;rowData`historicalPnl),enlist[groupColumn]!enlist groupVal;
        rIdx+:1];
    .var.riskContribution[pnlTable;groupColumn;confidenceLevel]
 };

.replay.replaySummary:{[replayResult]
    statusCol:replayResult`status;
    okCount:sum statusCol=`OK;
    errCount:(count replayResult)-okCount;
    pnlDist:.replay.historicalPnlDistribution replayResult;
    pnlVec:pnlDist`pnl;
    `scenarioCount`tradeRowCount`okRows`errorRows`minPnl`maxPnl`meanPnl`status!(
        count pnlDist;count replayResult;okCount;errCount;
        min pnlVec;max pnlVec;avg pnlVec;`OK)
 };
