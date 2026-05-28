/ timing.q - timing and performance utilities

.timing.timeFunction:{[functionHandle;argumentList]
    startTime:.z.p;
    functionResult:.[functionHandle;argumentList;{x}];
    endTime:.z.p;
    elapsedNanos:`long$endTime-startTime;
    elapsedMs:elapsedNanos%1000000;
    `result`elapsedMs`startTime`endTime!(functionResult;elapsedMs;startTime;endTime)
 };

.timing.benchmarkPricingGrid:{[tradeRow;marketData;pricingModel;configList]
    numConfigs:count configList;
    resultList:();
    loopIdx:0;
    while[loopIdx<numConfigs;
        currentConfig:configList loopIdx;
        timedResult:.timing.timeFunction[.engine.priceOption;(tradeRow;marketData;pricingModel;currentConfig)];
        pricedOk:99h=type timedResult`result;
        unitPrice:0Nf;
        if[pricedOk; unitPrice:(timedResult`result)`unitPrice];
        benchmarkRow:`method`numberOfSpotSteps`numberOfTimeSteps`elapsedMs`unitPrice`pricedOk!(
            currentConfig`method;
            currentConfig`numberOfSpotSteps;
            currentConfig`numberOfTimeSteps;
            timedResult`elapsedMs;
            unitPrice;
            pricedOk);
        resultList:resultList,enlist benchmarkRow;
        loopIdx+:1];
    resultList
 };

.timing.benchmarkPortfolioSize:{[tradeTable;marketDataBook;pricingModel;fdmConfig]
    portfolioSizes:1 2 4 8 16;
    resultList:();
    loopIdx:0;
    while[loopIdx<count portfolioSizes;
        portfolioSize:portfolioSizes loopIdx;
        subsetSize:portfolioSize&count tradeTable;
        subsetTable:subsetSize#tradeTable;
        timedResult:.timing.timeFunction[.portfolio.priceTradesWithMarketDataBook;(subsetTable;marketDataBook;pricingModel;fdmConfig)];
        benchmarkRow:`portfolioSize`elapsedMs!(subsetSize;timedResult`elapsedMs);
        resultList:resultList,enlist benchmarkRow;
        loopIdx+:1];
    resultList
 };
