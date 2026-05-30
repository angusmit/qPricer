/ electricity.q - electricity price foundation (v0.32)
/ No option pricing yet - only price replay, seasonality, spikes, spark spread

.commodity.electricity.validateHubPriceTable:{[hubPriceTable]
    if[0=count hubPriceTable; '"Empty hub price table"];
    firstRow:hubPriceTable 0;
    requiredCols:`snapshotDate`hub`price;
    tableCols:key firstRow;
    missingCols:requiredCols where not requiredCols in tableCols;
    if[0<count missingCols; '"Missing hub price columns: ",raze ", " ,/: string missingCols];
 };

.commodity.electricity.seasonalityFeatures:{[hubPriceTable]
    .commodity.electricity.validateHubPriceTable hubPriceTable;
    nRows:count hubPriceTable;
    dates:hubPriceTable`snapshotDate;
    monthVals:nRows#0i;
    dowVals:nRows#0i;
    sIdx:0;
    while[sIdx<nRows;
        monthVals[sIdx]:`mm$dates sIdx;
        dowVals[sIdx]:dates[sIdx] mod 7;
        sIdx+:1];
    hubPriceTable,'flip `monthOfYear`dayOfWeek!(monthVals;dowVals)
 };

.commodity.electricity.spikeFlags:{[hubPriceTable;thresholdConfig]
    prices:hubPriceTable`price;
    nRows:count prices;
    meanPx:avg prices;
    stdPx:dev prices;
    nStd:$[`nStd in key thresholdConfig;thresholdConfig`nStd;3f];
    spikeThreshold:meanPx+nStd*stdPx;
    spikeVals:prices>spikeThreshold;
    hubPriceTable,'flip `spikeFlag`spikeThreshold!(spikeVals;nRows#spikeThreshold)
 };

.commodity.electricity.monthlyAveragePrice:{[hubPriceTable]
    dates:hubPriceTable`snapshotDate;
    prices:hubPriceTable`price;
    hubs:hubPriceTable`hub;
    nRows:count hubPriceTable;
    monthKeys:nRows#0i;
    mkIdx:0;
    while[mkIdx<nRows;
        monthKeys[mkIdx]:12*(`year$dates mkIdx)+`mm$dates mkIdx;
        mkIdx+:1];
    uniqueMonths:asc distinct monthKeys;
    resultRows:();
    umIdx:0;
    while[umIdx<count uniqueMonths;
        mk:uniqueMonths umIdx;
        monthMask:monthKeys=mk;
        monthPrices:prices where monthMask;
        monthHub:first hubs where monthMask;
        monthDate:first dates where monthMask;
        resultRows:resultRows,enlist `monthKey`hub`snapshotDate`avgPrice`minPrice`maxPrice`observationCount!(
            mk;monthHub;monthDate;avg monthPrices;min monthPrices;max monthPrices;count monthPrices);
        umIdx+:1];
    resultRows
 };

.commodity.electricity.peakOffpeakProxy:{[hubPriceTable]
    / Simple proxy: weekday = peak (dow 1-5), weekend = off-peak (dow 0,6)
    nRows:count hubPriceTable;
    dates:hubPriceTable`snapshotDate;
    peakFlags:nRows#0b;
    poIdx:0;
    while[poIdx<nRows;
        dowVal:dates[poIdx] mod 7;
        peakFlags[poIdx]:dowVal within 1 5;
        poIdx+:1];
    hubPriceTable,'flip enlist[`peakFlag]!enlist peakFlags
 };

.commodity.electricity.sparkSpread:{[powerPx;gasPx;heatRate]
    powerPx-heatRate*gasPx
 };

.commodity.electricity.sparkSpreadTable:{[powerTable;gasTable;heatRate]
    .commodity.electricity.validateHubPriceTable powerTable;
    .commodity.electricity.validateHubPriceTable gasTable;
    nRows:count powerTable;
    sparkSpreads:nRows#0Nf;
    ssIdx:0;
    while[ssIdx<nRows;
        powerPx:(powerTable`price) ssIdx;
        / Find matching gas price by date
        gasDateMask:(gasTable`snapshotDate)=(powerTable`snapshotDate) ssIdx;
        gasRows:gasTable where gasDateMask;
        gasPx:$[0<count gasRows;(gasRows 0)`price;0Nf];
        sparkSpreads[ssIdx]:$[not null gasPx;powerPx-heatRate*gasPx;0Nf];
        ssIdx+:1];
    powerTable,'flip enlist[`sparkSpread]!enlist sparkSpreads
 };
