/ commodityFutures.q - futures/forward curve model (v0.32)

.commodity.futures.validateFuturesRow:{[futuresRow]
    if[not `settlementPrice in key futuresRow; '"Missing settlementPrice"];
    if[futuresRow[`settlementPrice]<=0f; '"Non-positive settlementPrice"];
    if[not `contractCode in key futuresRow; '"Missing contractCode"];
    if[not `commodity in key futuresRow; '"Missing commodity"];
 };

.commodity.futures.validateFuturesTable:{[futuresTable]
    if[0=count futuresTable; '"Empty futures table"];
    firstRow:futuresTable 0;
    requiredCols:`snapshotDate`commodity`contractCode`deliveryMonth`expiryDate`settlementPrice;
    tableCols:key firstRow;
    missingCols:requiredCols where not requiredCols in tableCols;
    if[0<count missingCols; '"Missing futures columns: ",raze ", " ,/: string missingCols];
 };

.commodity.futures.contractMonthFromCode:{[contractCode]
    codeStr:string contractCode;
    monthPart:(-2)#codeStr;
    "M"$monthPart
 };

.commodity.futures.buildCurveByDate:{[futuresTable]
    .commodity.futures.validateFuturesTable futuresTable;
    snapshotDates:distinct futuresTable`snapshotDate;
    resultRows:();
    sdIdx:0;
    while[sdIdx<count snapshotDates;
        sd:snapshotDates sdIdx;
        dayMask:(futuresTable`snapshotDate)=sd;
        dayRows:futuresTable where dayMask;
        / Sort by delivery month ascending
        sortIdx:iasc dayRows`deliveryMonth;
        nContracts:count sortIdx;
        tIdx:0;
        while[tIdx<nContracts;
            origIdx:sortIdx tIdx;
            rowData:dayRows origIdx;
            dte:rowData[`expiryDate]-sd;
            mltVal:$[`contractMultiplier in key rowData;rowData`contractMultiplier;1000f];
            statusVal:$[dte>=0i;`OK;`error];
            errMsg:$[dte<0i;"expired contract";""];
            resultRows:resultRows,enlist `snapshotDate`commodity`exchange`contractCode`deliveryMonth`expiryDate`tenorRank`settlementPrice`volume`openInterest`daysToExpiry`contractMultiplier`status`errorMessage!(
                sd;rowData`commodity;
                $[`exchange in key rowData;rowData`exchange;`NYMEX];
                rowData`contractCode;rowData`deliveryMonth;rowData`expiryDate;
                1+tIdx;rowData`settlementPrice;
                $[`volume in key rowData;rowData`volume;0f];
                $[`openInterest in key rowData;rowData`openInterest;0f];
                dte;mltVal;statusVal;errMsg);
            tIdx+:1];
        sdIdx+:1];
    resultRows
 };

.commodity.futures.frontContract:{[futuresCurve;snapshotDate]
    .commodity.futures.nthContract[futuresCurve;snapshotDate;1]
 };

.commodity.futures.nthContract:{[futuresCurve;snapshotDate;tenorN]
    dateMask:(futuresCurve`snapshotDate)=snapshotDate;
    okMask:(futuresCurve`status)=`OK;
    dayRows:futuresCurve where dateMask and okMask;
    if[0=count dayRows; '"No contracts on date ",string snapshotDate];
    tenorMask:(dayRows`tenorRank)=tenorN;
    matchRows:dayRows where tenorMask;
    if[0=count matchRows; '"No contract at tenor rank ",string tenorN];
    matchRows 0
 };

.commodity.futures.curveSlope:{[futuresCurve;snapshotDate;nearTenor;farTenor]
    nearRow:.commodity.futures.nthContract[futuresCurve;snapshotDate;nearTenor];
    farRow:.commodity.futures.nthContract[futuresCurve;snapshotDate;farTenor];
    farRow[`settlementPrice]-nearRow`settlementPrice
 };

.commodity.futures.rollSchedule:{[futuresCurve;rollDaysBefore]
    okMask:(futuresCurve`status)=`OK;
    frontMask:(futuresCurve`tenorRank)=1;
    frontRows:futuresCurve where okMask and frontMask;
    rollMask:(frontRows`daysToExpiry)<=rollDaysBefore;
    frontRows where rollMask
 };

.commodity.futures.futuresPnl:{[entryPx;exitPx;contractMult;qty]
    qty*contractMult*exitPx-entryPx
 };
