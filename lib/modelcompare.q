/ modelcompare.q - model comparison framework (v0.22)

.modelcompare.priceWithModel:{[modelName;optionTable;marketDataBook;calibrationResult;configDict]
    resultRows:();
    loopIdx:0;
    while[loopIdx<count optionTable;
        optionRow:optionTable loopIdx;
        rowResult:.modelcompare.__priceSingleOption[modelName;optionRow;marketDataBook;calibrationResult;configDict];
        resultRows:resultRows,enlist rowResult;
        loopIdx+:1];
    resultRows
 };

.modelcompare.__priceSingleOption:{[modelName;optionRow;marketDataBook;calibrationResult;configDict]
    errorRow:`modelName`optionId`underlying`optionType`strike`expiry`marketPrice`modelPrice`absoluteError`relativeError`squaredError`status`errorMessage!(
        modelName;optionRow`optionId;optionRow`underlying;optionRow`optionType;optionRow`strike;optionRow`expiry;optionRow`marketPrice;
        0Nf;0Nf;0Nf;0Nf;`ERROR;"");
    / BS surface: look up calibrated IV for matching option, compute BS price
    if[modelName~`blackScholesSurface;
        matchIdx:.modelcompare.__findOptionIdx[calibrationResult;optionRow`optionId];
        if[null matchIdx; :@[errorRow;`errorMessage;:;"No matching calibration row"]];
        calibRow:calibrationResult matchIdx;
        if[not calibRow[`status]~`OK; :@[errorRow;`errorMessage;:;"Calibration row failed"]];
        modelPriceVal:calibRow`modelPrice;
        :.modelcompare.__buildSuccessRow[modelName;optionRow;modelPriceVal]];
    / Heston grid: price with best Heston params
    if[modelName~`hestonGrid;
        bestParams:@[.calibration.bestCalibrationResult;calibrationResult;{x}];
        if[10h=type bestParams; :@[errorRow;`errorMessage;:;bestParams]];
        spotVal:.marketbook.getSpot[marketDataBook;optionRow`underlying];
        riskFreeRate:.marketbook.getRiskFreeRate[marketDataBook;optionRow`expiry];
        divYield:.marketbook.getDividendYield[marketDataBook;optionRow`underlying;optionRow`expiry];
        mktData:`underlying`spot`riskFreeRate`dividendYield`volatility!(optionRow`underlying;spotVal;riskFreeRate;divYield;0.2);
        hestonParams:`initialVariance`longRunVariance`meanReversion`volOfVol`correlation!(
            bestParams`initialVariance;bestParams`longRunVariance;bestParams`meanReversion;bestParams`volOfVol;bestParams`correlation);
        modelPriceVal:.calibration.priceWithHestonParams[optionRow;mktData;hestonParams;configDict];
        if[null modelPriceVal; :@[errorRow;`errorMessage;:;"Heston pricing failed"]];
        :.modelcompare.__buildSuccessRow[modelName;optionRow;modelPriceVal]];
    / SABR surface: price with SABR IV + BS
    if[modelName~`sabrSurface;
        spotVal:.marketbook.getSpot[marketDataBook;optionRow`underlying];
        riskFreeRate:.marketbook.getRiskFreeRate[marketDataBook;optionRow`expiry];
        divYield:.marketbook.getDividendYield[marketDataBook;optionRow`underlying;optionRow`expiry];
        mktData:`underlying`spot`riskFreeRate`dividendYield`volatility!(optionRow`underlying;spotVal;riskFreeRate;divYield;0.2);
        sabrResult:.sabr.priceEuropeanWithSabr[optionRow;mktData;calibrationResult;configDict];
        if[not sabrResult[`status]~`OK; :@[errorRow;`errorMessage;:;sabrResult`statusMessage]];
        :.modelcompare.__buildSuccessRow[modelName;optionRow;sabrResult`unitPrice]];
    / Merton series: price with Merton series formula
    if[modelName~`mertonSeries;
        spotVal:.marketbook.getSpot[marketDataBook;optionRow`underlying];
        riskFreeRate:.marketbook.getRiskFreeRate[marketDataBook;optionRow`expiry];
        divYield:.marketbook.getDividendYield[marketDataBook;optionRow`underlying;optionRow`expiry];
        mertonParams:calibrationResult,`riskFreeRate`dividendYield!(riskFreeRate;divYield);
        termCountValue:$[`termCount in key configDict;configDict`termCount;30];
        mertonPrice:@[.merton.priceEuropeanSeries[optionRow`optionType;spotVal;optionRow`strike;optionRow`expiry;mertonParams;];termCountValue;{x}];
        if[10h=type mertonPrice; :@[errorRow;`errorMessage;:;mertonPrice]];
        :.modelcompare.__buildSuccessRow[modelName;optionRow;mertonPrice]];
    @[errorRow;`errorMessage;:;"Unknown model: ",string modelName]
 };

.modelcompare.__findOptionIdx:{[calibrationResult;targetOptionId]
    optionIds:calibrationResult`optionId;
    matchIdx:optionIds?targetOptionId;
    if[matchIdx>=count optionIds; :0N];
    matchIdx
 };

.modelcompare.__buildSuccessRow:{[modelName;optionRow;modelPriceVal]
    mktPrice:optionRow`marketPrice;
    absErr:abs modelPriceVal-mktPrice;
    relErr:$[mktPrice>0f;absErr%mktPrice;0Nf];
    sqErr:absErr*absErr;
    `modelName`optionId`underlying`optionType`strike`expiry`marketPrice`modelPrice`absoluteError`relativeError`squaredError`status`errorMessage!(
        modelName;optionRow`optionId;optionRow`underlying;optionRow`optionType;optionRow`strike;optionRow`expiry;mktPrice;
        modelPriceVal;absErr;relErr;sqErr;`OK;"")
 };

.modelcompare.compareModels:{[modelList;optionTable;marketDataBook;calibrationDict;configDict]
    allRows:();
    modelIdx:0;
    while[modelIdx<count modelList;
        modelName:modelList modelIdx;
        calibResult:calibrationDict modelName;
        modelRows:.modelcompare.priceWithModel[modelName;optionTable;marketDataBook;calibResult;configDict];
        allRows:allRows,modelRows;
        modelIdx+:1];
    allRows
 };

.modelcompare.residualsByOption:{[comparisonTable]
    comparisonTable
 };

.modelcompare.rankModels:{[comparisonTable]
    modelNames:distinct comparisonTable`modelName;
    rankRows:();
    modelIdx:0;
    while[modelIdx<count modelNames;
        mn:modelNames modelIdx;
        modelMask:(comparisonTable`modelName)=mn;
        modelRows:comparisonTable where modelMask;
        okMask:(modelRows`status)=`OK;
        okRows:modelRows where okMask;
        pricedCount:count okRows;
        failedCount:(count modelRows)-pricedCount;
        rmseVal:0Nf; maeVal:0Nf; mreVal:0Nf; maxAeVal:0Nf; sseVal:0Nf;
        if[pricedCount>0;
            modelPrices:okRows`modelPrice;
            marketPrices:okRows`marketPrice;
            rmseVal:.objective.rmse[modelPrices;marketPrices];
            maeVal:.objective.mae[modelPrices;marketPrices];
            mreVal:.objective.meanRelativeError[modelPrices;marketPrices];
            maxAeVal:max okRows`absoluteError;
            diffs:modelPrices-marketPrices;
            sseVal:sum diffs*diffs];
        rankRows:rankRows,enlist `modelName`optionCount`pricedRows`failedRows`rmse`mae`meanRelativeError`maxAbsoluteError`totalSse`rank`status`errorMessage!(
            mn;count modelRows;pricedCount;failedCount;rmseVal;maeVal;mreVal;maxAeVal;sseVal;0N;`OK;"");
        modelIdx+:1];
    / Assign ranks by rmse
    rmseVals:rankRows`rmse;
    sortedIndices:iasc rmseVals;
    rankIdx:0;
    while[rankIdx<count sortedIndices;
        origIdx:sortedIndices rankIdx;
        rowDict:rankRows origIdx;
        rankRows[origIdx]:@[rowDict;`rank;:;1+rankIdx];
        rankIdx+:1];
    rankRows
 };

.modelcompare.bestModel:{[comparisonTable]
    ranked:.modelcompare.rankModels comparisonTable;
    ranks:{x`rank} each ranked;
    bestIdx:ranks?1;
    ranked bestIdx
 };

.modelcompare.applyBestModel:{[tradeTable;marketDataBook;comparisonTable;calibrationDict;configDict]
    bestModelRow:.modelcompare.bestModel comparisonTable;
    bestModelName:bestModelRow`modelName;
    calibResult:calibrationDict bestModelName;
    / For BS surface: extract a single calibration row for each trade (match by strike/expiry)
    / For Heston: use best params
    if[bestModelName~`blackScholesSurface;
        :.calibration.applyCalibratedModel[tradeTable;marketDataBook;calibResult 0;configDict]];
    if[bestModelName~`hestonGrid;
        bestParams:.calibration.bestCalibrationResult calibResult;
        :.calibration.applyCalibratedModel[tradeTable;marketDataBook;bestParams;configDict]];
    / Fallback
    resultList:();
    loopIdx:0;
    while[loopIdx<count tradeTable;
        resultList:resultList,enlist `tradeId`unitPrice`notionalPrice`status`statusMessage!(
            (tradeTable loopIdx)`tradeId;0Nf;0Nf;`ERROR;"Unsupported best model: ",string bestModelName);
        loopIdx+:1];
    resultList
 };
