/ calibration.q - model calibration framework (v0.21)

/ --- Validation ---

.calibration.validateMarketOptionTable:{[optionTable]
    rowCount:count optionTable;
    if[rowCount=0; '"Option table is empty"];
    requiredCols:`optionType`strike`expiry`marketPrice;
    firstRow:optionTable 0;
    tableCols:key firstRow;
    missingCols:requiredCols where not requiredCols in tableCols;
    if[0<count missingCols; '"Missing required columns: ",", " sv string missingCols];
    / Use column access (not each) since optionTable may be a proper table
    opTypes:optionTable`optionType;
    if[not all opTypes in `call`put; '"All optionType must be `call or `put"];
    strikes:optionTable`strike;
    if[not all strikes>0f; '"All strikes must be positive"];
    expiries:optionTable`expiry;
    if[not all expiries>0f; '"All expiries must be positive"];
    mktPrices:optionTable`marketPrice;
    if[not all mktPrices>0f; '"All marketPrice must be positive"];
 };

/ --- Black-Scholes IV calibration ---

.calibration.blackScholesVolForRow:{[optionRow;marketData;configDict]
    ivFn:{.iv.impliedVolatility[x 0;x 1;x 2;x 3;x 4;x 5;x 6]};
    ivArgs:(optionRow`optionType;optionRow`marketPrice;marketData`spot;optionRow`strike;optionRow`expiry;marketData`riskFreeRate;marketData`dividendYield);
    ivResult:@[ivFn;ivArgs;{x}];
    if[10h=type ivResult;
        :`optionId`underlying`optionType`strike`expiry`marketPrice`impliedVolatility`modelPrice`absoluteError`relativeError`status`errorMessage!(
            optionRow`optionId;optionRow`underlying;optionRow`optionType;optionRow`strike;optionRow`expiry;optionRow`marketPrice;
            0Nf;0Nf;0Nf;0Nf;`ERROR;ivResult)];
    calibratedIV:ivResult`impliedVolatility;
    modelPrice:.validation.blackScholesClosedForm[optionRow`optionType;marketData`spot;optionRow`strike;optionRow`expiry;marketData`riskFreeRate;marketData`dividendYield;calibratedIV];
    absErr:abs modelPrice-optionRow`marketPrice;
    relErr:$[optionRow[`marketPrice]>0f;absErr%optionRow`marketPrice;0Nf];
    `optionId`underlying`optionType`strike`expiry`marketPrice`impliedVolatility`modelPrice`absoluteError`relativeError`status`errorMessage!(
        optionRow`optionId;optionRow`underlying;optionRow`optionType;optionRow`strike;optionRow`expiry;optionRow`marketPrice;
        calibratedIV;modelPrice;absErr;relErr;`OK;"")
 };

.calibration.calibrateBlackScholesVolSurface:{[optionTable;marketDataBook;configDict]
    .calibration.validateMarketOptionTable optionTable;
    resultList:();
    loopIdx:0;
    while[loopIdx<count optionTable;
        optionRow:optionTable loopIdx;
        symbolName:optionRow`underlying;
        mktData:`underlying`spot`riskFreeRate`dividendYield`volatility!(
            symbolName;
            .marketbook.getSpot[marketDataBook;symbolName];
            .marketbook.getRiskFreeRate[marketDataBook;optionRow`expiry];
            .marketbook.getDividendYield[marketDataBook;symbolName;optionRow`expiry];
            0.2);
        rowResult:.calibration.blackScholesVolForRow[optionRow;mktData;configDict];
        resultList:resultList,enlist rowResult;
        loopIdx+:1];
    resultList
 };

/ --- Vol surface from calibration ---

.calibration.buildVolSurfaceFromCalibration:{[calibrationResult]
    statusCol:calibrationResult`status;
    okMask:statusCol=`OK;
    okRows:calibrationResult where okMask;
    surfaceRows:();
    loopIdx:0;
    while[loopIdx<count okRows;
        rowData:okRows loopIdx;
        surfaceRows:surfaceRows,enlist `underlying`strike`expiry`impliedVolatility!(
            rowData`underlying;rowData`strike;rowData`expiry;rowData`impliedVolatility);
        loopIdx+:1];
    surfaceRows
 };

/ --- Surface diagnostics ---

.calibration.surfaceDiagnostics:{[calibrationResult]
    statusCol:calibrationResult`status;
    okMask:statusCol=`OK;
    okRows:calibrationResult where okMask;
    failedCount:(count calibrationResult)-count okRows;
    if[0=count okRows;
        :`optionCount`averageIV`minIV`maxIV`averageAbsError`maxAbsError`failedRows!(
            0;0Nf;0Nf;0Nf;0Nf;0Nf;failedCount)];
    ivValues:okRows`impliedVolatility;
    absErrors:okRows`absoluteError;
    `optionCount`averageIV`minIV`maxIV`averageAbsError`maxAbsError`failedRows!(
        count okRows;avg ivValues;min ivValues;max ivValues;avg absErrors;max absErrors;failedCount)
 };

/ --- Heston parameter grid ---

.calibration.hestonParameterGrid:{[paramGridDict]
    ivList:paramGridDict`initialVarianceList;
    lvList:paramGridDict`longRunVarianceList;
    mrList:paramGridDict`meanReversionList;
    vvList:paramGridDict`volOfVolList;
    corrList:paramGridDict`correlationList;
    gridRows:();
    ivIdx:0;
    while[ivIdx<count ivList;
        lvIdx:0;
        while[lvIdx<count lvList;
            mrIdx:0;
            while[mrIdx<count mrList;
                vvIdx:0;
                while[vvIdx<count vvList;
                    corrIdx:0;
                    while[corrIdx<count corrList;
                        gridRows:gridRows,enlist `initialVariance`longRunVariance`meanReversion`volOfVol`correlation!(
                            ivList ivIdx;lvList lvIdx;mrList mrIdx;vvList vvIdx;corrList corrIdx);
                        corrIdx+:1];
                    vvIdx+:1];
                mrIdx+:1];
            lvIdx+:1];
        ivIdx+:1];
    gridRows
 };

/ --- Heston grid calibration ---

.calibration.priceWithHestonParams:{[optionRow;marketData;hestonParams;configDict]
    mcConfig:$[`mcConfig in key configDict; configDict`mcConfig; .montecarlo.defaultMcConfig[]];
    fullParams:hestonParams,`riskFreeRate`dividendYield!(marketData`riskFreeRate;marketData`dividendYield);
    trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
        optionRow`optionId;optionRow`underlying;`equityOption;`european;optionRow`optionType;optionRow`strike;optionRow`expiry;1f);
    hestonConfig:`mcConfig`hestonParams`modelType!(mcConfig;fullParams;`heston);
    hestonPriceFn:{.heston.priceEuropean[x 0;x 1;x 2]};
    priceResult:@[hestonPriceFn;(trade;marketData;hestonConfig);{x}];
    if[10h=type priceResult; :0Nf];
    priceResult`unitPrice
 };

.calibration.calibrateHestonGrid:{[optionTable;marketData;paramGridDict;configDict]
    .calibration.validateMarketOptionTable optionTable;
    paramGrid:.calibration.hestonParameterGrid paramGridDict;
    calibrationTable:();
    gridIdx:0;
    while[gridIdx<count paramGrid;
        paramRow:paramGrid gridIdx;
        modelPrices:();
        marketPrices:();
        failedPriceCount:0;
        optIdx:0;
        while[optIdx<count optionTable;
            optRow:optionTable optIdx;
            modelPriceVal:.calibration.priceWithHestonParams[optRow;marketData;paramRow;configDict];
            if[null modelPriceVal;
                failedPriceCount+:1];
            if[not null modelPriceVal;
                modelPrices:modelPrices,modelPriceVal;
                marketPrices:marketPrices,optRow`marketPrice];
            optIdx+:1];
        pricedCount:count modelPrices;
        rmseVal:$[pricedCount>0; .objective.rmse[modelPrices;marketPrices]; 0Nf];
        maeVal:$[pricedCount>0; .objective.mae[modelPrices;marketPrices]; 0Nf];
        sseVal:$[pricedCount>0; sum {x*x} each modelPrices-marketPrices; 0Nf];
        calibrationTable:calibrationTable,enlist `calibrationId`initialVariance`longRunVariance`meanReversion`volOfVol`correlation`rmse`mae`totalSse`pricedRows`failedRows`status`errorMessage!(
            gridIdx+1;paramRow`initialVariance;paramRow`longRunVariance;paramRow`meanReversion;paramRow`volOfVol;paramRow`correlation;
            rmseVal;maeVal;sseVal;pricedCount;failedPriceCount;`OK;"");
        gridIdx+:1];
    calibrationTable
 };

.calibration.bestCalibrationResult:{[calibrationTable]
    statusCol:calibrationTable`status;
    okMask:statusCol=`OK;
    okRows:calibrationTable where okMask;
    if[0=count okRows; '"No successful calibration results"];
    rmseValues:okRows`rmse;
    bestIdx:rmseValues?min rmseValues;
    okRows bestIdx
 };

/ --- Apply calibrated model ---

.calibration.applyCalibratedModel:{[tradeTable;marketDataBook;calibrationResult;configDict]
    resultList:();
    loopIdx:0;
    while[loopIdx<count tradeTable;
        tradeRow:tradeTable loopIdx;
        rowResult:.calibration.__priceWithCalibration[tradeRow;marketDataBook;calibrationResult;configDict];
        resultList:resultList,enlist rowResult;
        loopIdx+:1];
    resultList
 };

.calibration.__priceWithCalibration:{[tradeRow;marketDataBook;calibrationResult;configDict]
    / If calibrationResult has impliedVolatility field, it's BS surface calibration
    if[`impliedVolatility in key calibrationResult;
        iv:calibrationResult`impliedVolatility;
        spotVal:.marketbook.getSpot[marketDataBook;tradeRow`underlying];
        riskFreeRate:.marketbook.getRiskFreeRate[marketDataBook;tradeRow`expiry];
        divYield:.marketbook.getDividendYield[marketDataBook;tradeRow`underlying;tradeRow`expiry];
        bsFn:{.validation.blackScholesClosedForm[x 0;x 1;x 2;x 3;x 4;x 5;x 6]};
        modelPrice:@[bsFn;(tradeRow`optionType;spotVal;tradeRow`strike;tradeRow`expiry;riskFreeRate;divYield;iv);{x}];
        if[10h=type modelPrice;
            :`tradeId`unitPrice`notionalPrice`status`statusMessage!(tradeRow`tradeId;0Nf;0Nf;`ERROR;modelPrice)];
        :`tradeId`unitPrice`notionalPrice`status`statusMessage!(tradeRow`tradeId;modelPrice;modelPrice*tradeRow`notional;`OK;"")];
    / If calibrationResult has Heston params
    if[`initialVariance in key calibrationResult;
        hestonParams:`initialVariance`longRunVariance`meanReversion`volOfVol`correlation`riskFreeRate`dividendYield!(
            calibrationResult`initialVariance;calibrationResult`longRunVariance;calibrationResult`meanReversion;
            calibrationResult`volOfVol;calibrationResult`correlation;
            .marketbook.getRiskFreeRate[marketDataBook;tradeRow`expiry];
            .marketbook.getDividendYield[marketDataBook;tradeRow`underlying;tradeRow`expiry]);
        mcConfig:$[`mcConfig in key configDict; configDict`mcConfig; .montecarlo.defaultMcConfig[]];
        spotVal:.marketbook.getSpot[marketDataBook;tradeRow`underlying];
        mktData:`underlying`spot`riskFreeRate`dividendYield`volatility!(tradeRow`underlying;spotVal;hestonParams`riskFreeRate;hestonParams`dividendYield;0.2);
        hestonConfig:`mcConfig`hestonParams`modelType!(mcConfig;hestonParams;`heston);
        hestonFn:{.heston.priceEuropean[x 0;x 1;x 2]};
        priceResult:@[hestonFn;(tradeRow;mktData;hestonConfig);{x}];
        if[10h=type priceResult;
            :`tradeId`unitPrice`notionalPrice`status`statusMessage!(tradeRow`tradeId;0Nf;0Nf;`ERROR;priceResult)];
        :`tradeId`unitPrice`notionalPrice`status`statusMessage!(tradeRow`tradeId;priceResult`unitPrice;priceResult[`unitPrice]*tradeRow`notional;`OK;"")];
    `tradeId`unitPrice`notionalPrice`status`statusMessage!(tradeRow`tradeId;0Nf;0Nf;`ERROR;"Unknown calibration result format")
 };
