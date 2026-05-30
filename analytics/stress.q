/ stress.q - synthetic portfolio generation and stress testing

/ --- Public: data generation ---

.stress.generateSymbolList:{[symbolCountValue]
    `$"SYM" ,/: string 1+til symbolCountValue
 };

.stress.generateMarketDataBook:{[symbolCountValue;valuationDate]
    symbolList:.stress.generateSymbolList symbolCountValue;
    spotValues:50f+symbolCountValue?450f;
    volValues:0.15+symbolCountValue?0.30;
    divValues:symbolCountValue?0.02;
    spotTable:([] underlying:symbolList; spot:spotValues);
    volatilityTable:([] underlying:symbolList; volatility:volValues);
    rateTable:([] expiry:0.25 0.5 1 2f; riskFreeRate:0.03 0.04 0.05 0.055);
    dividendTable:([] underlying:symbolList; dividendYield:divValues);
    .marketbook.createMarketDataBook[spotTable;volatilityTable;rateTable;dividendTable]
 };

.stress.generateTradeTable:{[tradeCountValue;symbolList;valuationDate]
    tradeIds:1+til tradeCountValue;
    bookNames:tradeCountValue?`EQD`VOL`EXOTIC`HEDGE;
    underlyings:tradeCountValue?symbolList;
    optionTypes:tradeCountValue?`call`put;
    exerciseStyles:tradeCountValue?`european`european`european`european`american;
    strikes:50f+tradeCountValue?450f;
    expiries:0.25+tradeCountValue?1.75;
    notionals:(100000f+tradeCountValue?4900000f)*tradeCountValue?(1f;1f;1f;1f;-1f);
    ([] tradeId:tradeIds; bookName:bookNames; underlying:underlyings;
        productType:tradeCountValue#`equityOption;
        exerciseStyle:exerciseStyles; optionType:optionTypes;
        strike:strikes; expiry:expiries; notional:notionals;
        barrierType:tradeCountValue#`none;
        barrierLevel:tradeCountValue#0N;
        rebate:tradeCountValue#0f)
 };

.stress.generateSupportedTradeTable:{[tradeCountValue;symbolList;valuationDate]
    tradeIds:1+til tradeCountValue;
    bookNames:tradeCountValue?`EQD`VOL`EXOTIC`HEDGE;
    underlyings:tradeCountValue?symbolList;
    optionTypes:tradeCountValue?`call`put;
    / All European vanilla — guaranteed supported by CN and explicit
    strikes:50f+tradeCountValue?450f;
    expiries:0.25+tradeCountValue?1.75;
    notionals:(100000f+tradeCountValue?4900000f)*tradeCountValue?(1f;1f;1f;1f;-1f);
    ([] tradeId:tradeIds; bookName:bookNames; underlying:underlyings;
        productType:tradeCountValue#`equityOption;
        exerciseStyle:tradeCountValue#`european;
        optionType:optionTypes;
        strike:strikes; expiry:expiries; notional:notionals;
        barrierType:tradeCountValue#`none;
        barrierLevel:tradeCountValue#0N;
        rebate:tradeCountValue#0f)
 };

.stress.generateScenarioTable:{[]
    ([] scenario:`base`spotUp1Pct`spotDown1Pct`spotUp5Pct`spotDown5Pct`volUp1Point`volDown1Point`rateUp25Bp`rateDown25Bp;
        spotShiftPct:0 1 -1 5 -5 0 0 0 0f;
        volShiftAbs:0 0 0 0 0 0.01 -0.01 0 0f;
        rateShiftAbs:0 0 0 0 0 0 0 0.0025 -0.0025)
 };

.stress.generateLargeScenarioTable:{[]
    baseScenarios:.stress.generateScenarioTable[];
    extraScenarios:([]
        scenario:`spotUp10Pct`spotDown10Pct`spotUp20Pct`spotDown20Pct`volUp5Point`volDown5Point`rateUp100Bp`rateDown100Bp;
        spotShiftPct:10 -10 20 -20 0 0 0 0f;
        volShiftAbs:0 0 0 0 0.05 -0.05 0 0f;
        rateShiftAbs:0 0 0 0 0 0 0.01 -0.01);
    baseScenarios,extraScenarios
 };

/ --- Public: fault injection ---

.stress.injectMissingMarketData:{[marketDataBook;symbolsToRemove]
    sTable:marketDataBook`spotTable;
    vTable:marketDataBook`volatilityTable;
    dTable:marketDataBook`dividendTable;
    filteredSpot:sTable where not sTable[`underlying] in symbolsToRemove;
    filteredVol:vTable where not vTable[`underlying] in symbolsToRemove;
    filteredDiv:dTable where not dTable[`underlying] in symbolsToRemove;
    `spotTable`volatilityTable`rateTable`dividendTable!(filteredSpot;filteredVol;marketDataBook`rateTable;filteredDiv)
 };

.stress.injectBadTrades:{[tradeTable;badTradeCountValue]
    badIds:(count tradeTable)+1+til badTradeCountValue;
    badTable:([]
        tradeId:badIds;
        bookName:badTradeCountValue#`BAD;
        underlying:badTradeCountValue?`NOSYM1`NOSYM2;
        productType:badTradeCountValue#`equityOption;
        exerciseStyle:badTradeCountValue#`european;
        optionType:badTradeCountValue#`call;
        strike:badTradeCountValue#-100f;
        expiry:badTradeCountValue#0f;
        notional:badTradeCountValue#1000000f;
        barrierType:badTradeCountValue#`none;
        barrierLevel:badTradeCountValue#0N;
        rebate:badTradeCountValue#0f);
    tradeTable,badTable
 };

/ --- Public: stress runners ---

.stress.runPortfolioStress:{[tradeCountValue;symbolCountValue;configDict;modeSym]
    marketDataBook:.stress.generateMarketDataBook[symbolCountValue;2025.01.01];
    symbolList:marketDataBook[`spotTable]`underlying;
    tradeTable:$[modeSym~`supportedOnly;
        .stress.generateSupportedTradeTable[tradeCountValue;symbolList;2025.01.01];
        .stress.generateTradeTable[tradeCountValue;symbolList;2025.01.01]];
    previousBook:.stress.generateMarketDataBook[symbolCountValue;2024.12.31];
    startTime:.z.p;
    runResult:.batch.runDailyPricing[tradeTable;marketDataBook;previousBook;configDict];
    elapsedMs:(`long$.z.p-startTime)%1000000;
    pricingResult:runResult`pricingResult;
    okRows:sum pricingResult[`status]=`OK;
    errorRows:(count pricingResult)-okRows;
    `tradeCount`symbolCount`elapsedMs`okRows`errorRows`runResult!(
        tradeCountValue;symbolCountValue;elapsedMs;okRows;errorRows;runResult)
 };

.stress.runMissingDataStress:{[tradeCountValue;symbolCountValue;missingSymbolCount;configDict]
    fullBook:.stress.generateMarketDataBook[symbolCountValue;2025.01.01];
    symbolList:fullBook[`spotTable]`underlying;
    tradeTable:.stress.generateSupportedTradeTable[tradeCountValue;symbolList;2025.01.01];
    symbolsToRemove:missingSymbolCount#symbolList;
    damagedBook:.stress.injectMissingMarketData[fullBook;symbolsToRemove];
    pricingModel:configDict`model;
    fdmConfig:configDict`fdmConfig;
    pricingResult:.portfolio.priceTradesWithMarketDataBook[tradeTable;damagedBook;pricingModel;fdmConfig];
    okRows:sum pricingResult[`status]=`OK;
    errorRows:(count pricingResult)-okRows;
    `tradeCount`missingSymbols`okRows`errorRows`pricingResult!(
        tradeCountValue;missingSymbolCount;okRows;errorRows;pricingResult)
 };

.stress.runPnlExplainStress:{[tradeCountValue;symbolCountValue;configDict]
    book0:.stress.generateMarketDataBook[symbolCountValue;2025.01.01];
    book1:.stress.generateMarketDataBook[symbolCountValue;2025.01.02];
    symbolList:book0[`spotTable]`underlying;
    tradeTable:.stress.generateSupportedTradeTable[tradeCountValue;symbolList;2025.01.01];
    pnlConfig:`model`fdmConfig`timeStepYears`bookName!(configDict`model;configDict`fdmConfig;1%252;"stressBook");
    startTime:.z.p;
    explainResult:.pnl.explainPortfolio[tradeTable;book0;book1;pnlConfig];
    elapsedMs:(`long$.z.p-startTime)%1000000;
    okRows:sum explainResult[`status]=`OK;
    errorRows:(count explainResult)-okRows;
    `tradeCount`elapsedMs`okRows`errorRows`explainResult!(
        count tradeTable;elapsedMs;okRows;errorRows;explainResult)
 };

/ --- v0.14 extended generators ---

.stress.generateExtendedSupportedTradeTable:{[tradeCountValue;symbolList;valuationDate]
    / Mix: European vanilla, American call/put, European barrier (all 8)
    tradeIds:1+til tradeCountValue;
    bookNames:tradeCountValue?`EQD`VOL`EXOTIC`HEDGE;
    underlyings:tradeCountValue?symbolList;
    / 40% European vanilla, 30% American, 30% European barrier
    exercisePool:`european`european`european`european`american`american`american`european`european`european;
    barrierPool:`none`none`none`none`none`none`none`upAndOut`downAndOut`upAndIn;
    poolIndices:tradeCountValue?til count exercisePool;
    exerciseStyles:exercisePool poolIndices;
    barrierTypes:barrierPool poolIndices;
    optionTypes:tradeCountValue?`call`put;
    strikes:50f+tradeCountValue?450f;
    expiries:0.25+tradeCountValue?1.75;
    notionals:(100000f+tradeCountValue?4900000f)*tradeCountValue?(1f;1f;1f;1f;-1f);
    / Barrier levels: up=strike*1.3, down=strike*0.7
    barrierLevels:tradeCountValue#0N;
    barrierIdx:0;
    while[barrierIdx<tradeCountValue;
        barrierSym:barrierTypes barrierIdx;
        if[(barrierSym~`upAndOut) or barrierSym~`upAndIn;
            barrierLevels[barrierIdx]:`long$strikes[barrierIdx]*1.3];
        if[(barrierSym~`downAndOut) or barrierSym~`downAndIn;
            barrierLevels[barrierIdx]:`long$strikes[barrierIdx]*0.7];
        barrierIdx+:1];
    ([] tradeId:tradeIds; bookName:bookNames; underlying:underlyings;
        productType:tradeCountValue#`equityOption;
        exerciseStyle:exerciseStyles; optionType:optionTypes;
        strike:strikes; expiry:expiries; notional:notionals;
        barrierType:barrierTypes;
        barrierLevel:barrierLevels;
        rebate:tradeCountValue#0f)
 };

.stress.generateBarrierTradeTable:{[tradeCountValue;symbolList;valuationDate]
    tradeIds:1+til tradeCountValue;
    underlyings:tradeCountValue?symbolList;
    optionTypes:tradeCountValue?`call`put;
    barrierTypes:tradeCountValue?`upAndOut`downAndOut`upAndIn`downAndIn;
    strikes:50f+tradeCountValue?450f;
    expiries:0.25+tradeCountValue?1.75;
    notionals:100000f+tradeCountValue?4900000f;
    barrierLevels:tradeCountValue#0N;
    barrierIdx:0;
    while[barrierIdx<tradeCountValue;
        bSym:barrierTypes barrierIdx;
        if[(bSym~`upAndOut) or bSym~`upAndIn;
            barrierLevels[barrierIdx]:`long$strikes[barrierIdx]*1.3];
        if[(bSym~`downAndOut) or bSym~`downAndIn;
            barrierLevels[barrierIdx]:`long$strikes[barrierIdx]*0.7];
        barrierIdx+:1];
    ([] tradeId:tradeIds; bookName:tradeCountValue#`EXOTIC; underlying:underlyings;
        productType:tradeCountValue#`equityOption;
        exerciseStyle:tradeCountValue#`european;
        optionType:optionTypes;
        strike:strikes; expiry:expiries; notional:notionals;
        barrierType:barrierTypes; barrierLevel:barrierLevels; rebate:tradeCountValue#0f)
 };

.stress.generateAmericanTradeTable:{[tradeCountValue;symbolList;valuationDate]
    tradeIds:1+til tradeCountValue;
    underlyings:tradeCountValue?symbolList;
    optionTypes:tradeCountValue?`call`put;
    strikes:50f+tradeCountValue?450f;
    expiries:0.25+tradeCountValue?1.75;
    notionals:100000f+tradeCountValue?4900000f;
    ([] tradeId:tradeIds; bookName:tradeCountValue#`EQD; underlying:underlyings;
        productType:tradeCountValue#`equityOption;
        exerciseStyle:tradeCountValue#`american;
        optionType:optionTypes;
        strike:strikes; expiry:expiries; notional:notionals;
        barrierType:tradeCountValue#`none; barrierLevel:tradeCountValue#0N; rebate:tradeCountValue#0f)
 };

.stress.generateUnsupportedTradeTable:{[tradeCountValue]
    tradeIds:1+til tradeCountValue;
    / Intentionally invalid: American + barrier (rejected by product validation)
    ([] tradeId:tradeIds; bookName:tradeCountValue#`BAD; underlying:tradeCountValue#`NOSYM;
        productType:tradeCountValue#`equityOption;
        exerciseStyle:tradeCountValue#`american;
        optionType:tradeCountValue#`call;
        strike:tradeCountValue#100f; expiry:tradeCountValue#1f; notional:tradeCountValue#1000000f;
        barrierType:tradeCountValue#`upAndOut; barrierLevel:tradeCountValue#130; rebate:tradeCountValue#0f)
 };
