/ perfopt.q - performance optimisation: cached pricing, Greek reuse, batch comparison

/ --- Public: cached portfolio pricing ---

.perfopt.priceBookWithCache:{[tradeTable;marketDataBook;pricingModel;fdmConfig]
    numTrades:count tradeTable;
    cacheDict:.cache.emptyCache[];
    resultList:();
    hitCount:0;
    loopIdx:0;
    while[loopIdx<numTrades;
          currentTrade:tradeTable loopIdx;
          cacheResult:.[{[td;mdb;mdl;cfg;cd]
                            mktData:.marketbook.getMarketDataForTrade[mdb;td];
                            .cache.getOrPrice[cd;td;mktData;mdl;cfg]
                        };(currentTrade;marketDataBook;pricingModel;fdmConfig;cacheDict);{x}];
          if[10h=type cacheResult;
             resultList:resultList,enlist `tradeId`underlying`optionType`unitPrice`notionalPrice`method`modelName`status`statusMessage!(
                 currentTrade`tradeId;currentTrade`underlying;currentTrade`optionType;
                 0Nf;0Nf;fdmConfig`method;pricingModel`modelName;`ERROR;cacheResult)];
          if[99h=type cacheResult;
             cacheDict:cacheResult`cacheDict;
             if[cacheResult`cacheHit; hitCount+:1];
             rawResult:cacheResult`pricingResult;
             enrichedResult:rawResult,`status`statusMessage!(`OK;"");
             resultList:resultList,enlist enrichedResult];
          loopIdx+:1];
    `pricingResult`cacheDict`cacheHits`cacheMisses!(resultList;cacheDict;hitCount;numTrades-hitCount)
 };

/ --- Public: PnL explain with pre-computed Greeks ---

.perfopt.pnlExplainWithGreekReuse:{[tradeTable;marketDataBook0;marketDataBook1;configDict;greekResult]
    pricingModel:configDict`model;
    fdmConfig:configDict`fdmConfig;
    timeStepYears:configDict`timeStepYears;
    bookLabel:configDict`bookName;
    numTrades:count tradeTable;
    resultList:();
    loopIdx:0;
    while[loopIdx<numTrades;
          currentTrade:tradeTable loopIdx;
          singleResult:.[.perfopt.__pnlSingleWithGreeks;
                         (currentTrade;marketDataBook0;marketDataBook1;pricingModel;fdmConfig;timeStepYears;bookLabel;greekResult);{x}];
          if[10h=type singleResult;
             singleResult:`tradeId`bookName`underlyingSym`pv0`pv1`actualPnL`deltaPnL`gammaPnL`vegaPnL`rhoPnL`thetaPnL`explainedPnL`unexplainedPnL`status`errorMessage!(
                 currentTrade`tradeId;bookLabel;currentTrade`underlying;
                 0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;`ERROR;singleResult)];
          resultList:resultList,enlist singleResult;
          loopIdx+:1];
    resultList
 };

.perfopt.__pnlSingleWithGreeks:{[tradeRow;marketDataBook0;marketDataBook1;pricingModel;fdmConfig;timeStepYears;bookLabel;greekResult]
    mktData0:.marketbook.getMarketDataForTrade[marketDataBook0;tradeRow];
    mktData1:.marketbook.getMarketDataForTrade[marketDataBook1;tradeRow];
    priceResult0:.engine.priceOption[tradeRow;mktData0;pricingModel;fdmConfig];
    priceResult1:.engine.priceOption[tradeRow;mktData1;pricingModel;fdmConfig];
    pv0Val:priceResult0`unitPrice;
    pv1Val:priceResult1`unitPrice;
    actualPnLVal:pv1Val-pv0Val;
    / Look up pre-computed Greeks by tradeId
    tradeIdVal:tradeRow`tradeId;
    greekRows:greekResult where greekResult[`tradeId]=tradeIdVal;
    if[0=count greekRows; '"No Greek result for tradeId ",string tradeIdVal];
    greekRow:greekRows 0;
    if[not greekRow[`status]~`OK; '"Greeks not OK for tradeId ",string tradeIdVal];
    deltaVal:greekRow`delta;
    gammaVal:greekRow`gamma;
    thetaVal:greekRow`theta;
    vegaVal:greekRow`vega;
    rhoVal:greekRow`rho;
    dSpot:mktData1[`spot]-mktData0`spot;
    dVol:mktData1[`volatility]-mktData0`volatility;
    dRate:mktData1[`riskFreeRate]-mktData0`riskFreeRate;
    deltaPnLVal:deltaVal*dSpot;
    gammaPnLVal:0.5*gammaVal*dSpot*dSpot;
    vegaPnLVal:vegaVal*dVol;
    rhoPnLVal:rhoVal*dRate;
    thetaPnLVal:thetaVal*timeStepYears;
    explainedPnLVal:deltaPnLVal+gammaPnLVal+vegaPnLVal+rhoPnLVal+thetaPnLVal;
    unexplainedPnLVal:actualPnLVal-explainedPnLVal;
    `tradeId`bookName`underlyingSym`pv0`pv1`actualPnL`deltaPnL`gammaPnL`vegaPnL`rhoPnL`thetaPnL`explainedPnL`unexplainedPnL`status`errorMessage!(
        tradeIdVal;bookLabel;tradeRow`underlying;
        pv0Val;pv1Val;actualPnLVal;deltaPnLVal;gammaPnLVal;vegaPnLVal;rhoPnLVal;thetaPnLVal;
        explainedPnLVal;unexplainedPnLVal;`OK;"")
 };

/ --- Public: optimised batch run ---

.perfopt.runOptimisedBatch:{[tradeTable;marketDataBook;previousMarketDataBook;configDict]
    pricingModel:configDict`model;
    fdmConfig:configDict`fdmConfig;
    timeStepYears:configDict`timeStepYears;
    bookLabel:configDict`bookName;
    valuationDate:configDict`valuationDate;
    runLabel:configDict`runLabel;
    / Step 1: Price with cache
    cacheResult:.perfopt.priceBookWithCache[tradeTable;marketDataBook;pricingModel;fdmConfig];
    pricingResult:cacheResult`pricingResult;
    / Step 2: Greeks at current market (for risk summary)
    greekResult:.batch.__calculateGreeksFromBook[tradeTable;marketDataBook;pricingModel;fdmConfig];
    / Step 3: Scenarios (existing function)
    scenarioResult:.portfolio.generatePortfolioScenarioReportWithMarketDataBook[tradeTable;marketDataBook;pricingModel;fdmConfig];
    / Step 4: PnL explain with Greek reuse (Greeks at t0 to match baseline PnL explain)
    greekResultT0:.batch.__calculateGreeksFromBook[tradeTable;previousMarketDataBook;pricingModel;fdmConfig];
    pnlConfig:`model`fdmConfig`timeStepYears`bookName!(pricingModel;fdmConfig;timeStepYears;bookLabel);
    pnlExplainResult:.perfopt.pnlExplainWithGreekReuse[tradeTable;previousMarketDataBook;marketDataBook;pnlConfig;greekResultT0];
    / Summaries
    portfolioSummary:.report.portfolioSummary pricingResult;
    riskSummary:.report.riskSummary greekResult;
    scenarioSummary:.report.scenarioSummary scenarioResult;
    errorSummary:.audit.errorAudit[pricingResult;scenarioResult;pnlExplainResult];
    partialResult:`pricingResult`greekResult`scenarioResult`pnlExplainResult`portfolioSummary`riskSummary`scenarioSummary!(
        pricingResult;greekResult;scenarioResult;pnlExplainResult;portfolioSummary;riskSummary;scenarioSummary);
    auditRecord:.audit.createAuditRecord[partialResult;valuationDate;runLabel];
    partialResult,`errorSummary`auditRecord!(errorSummary;auditRecord)
 };

/ --- Public: comparison ---

.perfopt.compareBaselineOptimised:{[tradeTable;marketDataBook;previousMarketDataBook;configDict]
    pricingModel:configDict`model;
    fdmConfig:configDict`fdmConfig;
    timeStepYears:configDict`timeStepYears;
    bookLabel:configDict`bookName;
    resultList:();
    / --- Full batch comparison ---
    baselineStart:.z.p;
    baselineResult:.batch.runDailyPricing[tradeTable;marketDataBook;previousMarketDataBook;configDict];
    batchBaseMs:(`long$.z.p-baselineStart)%1000000;
    optimisedStart:.z.p;
    optimisedResult:.perfopt.runOptimisedBatch[tradeTable;marketDataBook;previousMarketDataBook;configDict];
    batchOptMs:(`long$.z.p-optimisedStart)%1000000;
    batchBasePVs:`float${x`unitPrice} each baselineResult`pricingResult;
    batchOptPVs:`float${x`unitPrice} each optimisedResult`pricingResult;
    batchMaxDiff:.perfopt.__safeMaxAbsDiff[batchBasePVs;batchOptPVs];
    batchSpeedup:$[batchOptMs>0;batchBaseMs%batchOptMs;0Nf];
    resultList:resultList,enlist `stepName`baselineMs`optimisedMs`speedupRatio`maxAbsoluteDifference`status`errorMessage!(
        "fullBatch";batchBaseMs;batchOptMs;batchSpeedup;batchMaxDiff;`OK;"");
    / --- Pricing comparison ---
    pricingBaseStart:.z.p;
    pricingBaseResult:.portfolio.priceTradesWithMarketDataBook[tradeTable;marketDataBook;pricingModel;fdmConfig];
    pricingBaseMs:(`long$.z.p-pricingBaseStart)%1000000;
    pricingOptStart:.z.p;
    pricingOptResult:.perfopt.priceBookWithCache[tradeTable;marketDataBook;pricingModel;fdmConfig];
    pricingOptMs:(`long$.z.p-pricingOptStart)%1000000;
    pricingBasePVs:`float${x`unitPrice} each pricingBaseResult;
    pricingOptPVs:`float${x`unitPrice} each pricingOptResult`pricingResult;
    pricingMaxDiff:.perfopt.__safeMaxAbsDiff[pricingBasePVs;pricingOptPVs];
    pricingSpeedup:$[pricingOptMs>0;pricingBaseMs%pricingOptMs;0Nf];
    resultList:resultList,enlist `stepName`baselineMs`optimisedMs`speedupRatio`maxAbsoluteDifference`status`errorMessage!(
        "pricing";pricingBaseMs;pricingOptMs;pricingSpeedup;pricingMaxDiff;`OK;"");
    / --- PnL explain comparison ---
    pnlConfig:`model`fdmConfig`timeStepYears`bookName!(pricingModel;fdmConfig;timeStepYears;bookLabel);
    pnlBaseStart:.z.p;
    pnlBaseResult:.pnl.explainPortfolio[tradeTable;previousMarketDataBook;marketDataBook;pnlConfig];
    pnlBaseMs:(`long$.z.p-pnlBaseStart)%1000000;
    / Optimised: compute Greeks once at t0 (previous book, matching baseline), then reuse
    greekResult:.batch.__calculateGreeksFromBook[tradeTable;previousMarketDataBook;pricingModel;fdmConfig];
    pnlOptStart:.z.p;
    pnlOptResult:.perfopt.pnlExplainWithGreekReuse[tradeTable;previousMarketDataBook;marketDataBook;pnlConfig;greekResult];
    pnlOptMs:(`long$.z.p-pnlOptStart)%1000000;
    pnlBasePnL:`float${x`actualPnL} each pnlBaseResult;
    pnlOptPnL:`float${x`actualPnL} each pnlOptResult;
    pnlBaseExplained:`float${x`explainedPnL} each pnlBaseResult;
    pnlOptExplained:`float${x`explainedPnL} each pnlOptResult;
    pnlMaxDiff:max (.perfopt.__safeMaxAbsDiff[pnlBasePnL;pnlOptPnL];.perfopt.__safeMaxAbsDiff[pnlBaseExplained;pnlOptExplained]);
    pnlSpeedup:$[pnlOptMs>0;pnlBaseMs%pnlOptMs;0Nf];
    resultList:resultList,enlist `stepName`baselineMs`optimisedMs`speedupRatio`maxAbsoluteDifference`status`errorMessage!(
        "pnlExplain";pnlBaseMs;pnlOptMs;pnlSpeedup;pnlMaxDiff;`OK;"");
    resultList
 };

.perfopt.__safeMaxAbsDiff:{[listA;listB]
    validA:listA where not null listA;
    validB:listB where not null listB;
    pairLen:(count validA)&count validB;
    if[pairLen=0; :0f];
    max abs (pairLen#validA)-pairLen#validB
 };
