/ calibreport.q - calibration reporting and bucket diagnostics (v0.22)

/ --- Residual report ---

.calibreport.optionResidualReport:{[comparisonTable]
    residuals:(comparisonTable`modelPrice)-comparisonTable`marketPrice;
    resultRows:();
    loopIdx:0;
    while[loopIdx<count comparisonTable;
        rowData:comparisonTable loopIdx;
        residualVal:$[rowData[`status]~`OK;rowData[`modelPrice]-rowData`marketPrice;0Nf];
        resultRows:resultRows,enlist `modelName`optionId`underlying`optionType`strike`expiry`marketPrice`modelPrice`residual`absoluteError`relativeError`status`errorMessage!(
            rowData`modelName;rowData`optionId;rowData`underlying;rowData`optionType;rowData`strike;rowData`expiry;
            rowData`marketPrice;rowData`modelPrice;residualVal;rowData`absoluteError;rowData`relativeError;
            rowData`status;rowData`errorMessage);
        loopIdx+:1];
    resultRows
 };

/ --- Maturity bucket report ---

.calibreport.__maturityBucket:{[expiryVal]
    if[expiryVal<=0.25; :`short];
    if[expiryVal<=1.0; :`medium];
    `long
 };

.calibreport.maturityBucketReport:{[comparisonTable]
    modelNames:distinct comparisonTable`modelName;
    bucketNames:`short`medium`long;
    resultRows:();
    modelIdx:0;
    while[modelIdx<count modelNames;
        mn:modelNames modelIdx;
        modelMask:(comparisonTable`modelName)=mn;
        modelRows:comparisonTable where modelMask;
        okMask:(modelRows`status)=`OK;
        okRows:modelRows where okMask;
        bucketIdx:0;
        while[bucketIdx<count bucketNames;
            bn:bucketNames bucketIdx;
            bucketMask:{.calibreport.__maturityBucket[x]~y}[;bn] each okRows`expiry;
            bucketRows:okRows where bucketMask;
            bucketCount:count bucketRows;
            rmseVal:0Nf; maeVal:0Nf; mreVal:0Nf; maxAeVal:0Nf;
            if[bucketCount>0;
                bModelPrices:bucketRows`modelPrice;
                bMarketPrices:bucketRows`marketPrice;
                rmseVal:.objective.rmse[bModelPrices;bMarketPrices];
                maeVal:.objective.mae[bModelPrices;bMarketPrices];
                mreVal:.objective.meanRelativeError[bModelPrices;bMarketPrices];
                maxAeVal:max bucketRows`absoluteError];
            resultRows:resultRows,enlist `modelName`maturityBucket`optionCount`rmse`mae`meanRelativeError`maxAbsoluteError!(
                mn;bn;bucketCount;rmseVal;maeVal;mreVal;maxAeVal);
            bucketIdx+:1];
        modelIdx+:1];
    resultRows
 };

/ --- Moneyness bucket report ---

.calibreport.__moneynessBucket:{[moneynessVal]
    if[moneynessVal<0.8; :`lowStrike];
    if[moneynessVal<0.95; :`belowAtm];
    if[moneynessVal<=1.05; :`atm];
    if[moneynessVal<=1.2; :`aboveAtm];
    `highStrike
 };

.calibreport.moneynessBucketReport:{[comparisonTable;marketDataBook]
    modelNames:distinct comparisonTable`modelName;
    bucketNames:`lowStrike`belowAtm`atm`aboveAtm`highStrike;
    resultRows:();
    modelIdx:0;
    while[modelIdx<count modelNames;
        mn:modelNames modelIdx;
        modelMask:(comparisonTable`modelName)=mn;
        modelRows:comparisonTable where modelMask;
        okMask:(modelRows`status)=`OK;
        okRows:modelRows where okMask;
        / Compute moneyness for each ok row
        moneynessVals:();
        loopIdx:0;
        while[loopIdx<count okRows;
            rowData:okRows loopIdx;
            spotVal:.marketbook.getSpot[marketDataBook;rowData`underlying];
            moneynessVals:moneynessVals,rowData[`strike]%spotVal;
            loopIdx+:1];
        bucketIdx:0;
        while[bucketIdx<count bucketNames;
            bn:bucketNames bucketIdx;
            bucketMask:{.calibreport.__moneynessBucket[x]~y}[;bn] each moneynessVals;
            bucketRows:okRows where bucketMask;
            bucketCount:count bucketRows;
            rmseVal:0Nf; maeVal:0Nf; mreVal:0Nf; maxAeVal:0Nf;
            if[bucketCount>0;
                bModelPrices:bucketRows`modelPrice;
                bMarketPrices:bucketRows`marketPrice;
                rmseVal:.objective.rmse[bModelPrices;bMarketPrices];
                maeVal:.objective.mae[bModelPrices;bMarketPrices];
                mreVal:.objective.meanRelativeError[bModelPrices;bMarketPrices];
                maxAeVal:max bucketRows`absoluteError];
            resultRows:resultRows,enlist `modelName`moneynessBucket`optionCount`rmse`mae`meanRelativeError`maxAbsoluteError!(
                mn;bn;bucketCount;rmseVal;maeVal;mreVal;maxAeVal);
            bucketIdx+:1];
        modelIdx+:1];
    resultRows
 };

/ --- Strike bucket report ---

.calibreport.strikeBucketReport:{[comparisonTable]
    modelNames:distinct comparisonTable`modelName;
    resultRows:();
    modelIdx:0;
    while[modelIdx<count modelNames;
        mn:modelNames modelIdx;
        modelMask:(comparisonTable`modelName)=mn;
        modelRows:comparisonTable where modelMask;
        okMask:(modelRows`status)=`OK;
        okRows:modelRows where okMask;
        if[0<count okRows;
            strikeVals:distinct okRows`strike;
            strikeIdx:0;
            while[strikeIdx<count strikeVals;
                strikeVal:strikeVals strikeIdx;
                strikeMask:(okRows`strike)=strikeVal;
                strikeRows:okRows where strikeMask;
                bModelPrices:strikeRows`modelPrice;
                bMarketPrices:strikeRows`marketPrice;
                resultRows:resultRows,enlist `modelName`strike`optionCount`rmse`mae`maxAbsoluteError!(
                    mn;strikeVal;count strikeRows;
                    .objective.rmse[bModelPrices;bMarketPrices];
                    .objective.mae[bModelPrices;bMarketPrices];
                    max strikeRows`absoluteError);
                strikeIdx+:1]];
        modelIdx+:1];
    resultRows
 };

/ --- Model summary ---

.calibreport.modelSummaryReport:{[comparisonTable]
    .modelcompare.rankModels comparisonTable
 };

/ --- CSV export ---

.calibreport.exportCalibrationReports:{[comparisonTable;outputDirectory;reportLabel]
    csvExportFn:{.report.exportCsv[x 0;x 1];`OK};
    exportResults:();
    / Residuals
    residualTable:.calibreport.optionResidualReport comparisonTable;
    residualPath:outputDirectory,"/",reportLabel,"_optionResiduals.csv";
    residualStatus:@[csvExportFn;(residualTable;residualPath);{`ERROR}];
    exportResults:exportResults,enlist `reportName`filePath`status`errorMessage!("optionResiduals";residualPath;residualStatus;"");
    / Model summary
    summaryTable:.calibreport.modelSummaryReport comparisonTable;
    summaryPath:outputDirectory,"/",reportLabel,"_modelSummary.csv";
    summaryStatus:@[csvExportFn;(summaryTable;summaryPath);{`ERROR}];
    exportResults:exportResults,enlist `reportName`filePath`status`errorMessage!("modelSummary";summaryPath;summaryStatus;"");
    / Maturity buckets
    matBuckets:.calibreport.maturityBucketReport comparisonTable;
    matPath:outputDirectory,"/",reportLabel,"_maturityBuckets.csv";
    matStatus:@[csvExportFn;(matBuckets;matPath);{`ERROR}];
    exportResults:exportResults,enlist `reportName`filePath`status`errorMessage!("maturityBuckets";matPath;matStatus;"");
    / Strike buckets
    strikeBuckets:.calibreport.strikeBucketReport comparisonTable;
    strikePath:outputDirectory,"/",reportLabel,"_strikeBuckets.csv";
    strikeStatus:@[csvExportFn;(strikeBuckets;strikePath);{`ERROR}];
    exportResults:exportResults,enlist `reportName`filePath`status`errorMessage!("strikeBuckets";strikePath;strikeStatus;"");
    exportResults
 };