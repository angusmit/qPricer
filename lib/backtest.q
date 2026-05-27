/ backtest.q — Barchart multi-day replay & performance analytics
/ Namespace: .parser.barchart (public) / .parser.barchart.__ (private)
/ Depends on: lib/parser.q (loaded via init.q)

/ ══════════════════════════════════════════════════════════════════
/ PUBLIC API
/ ══════════════════════════════════════════════════════════════════

/ Sorted unique snapshot dates where status=OK
.parser.barchart.availableDates:{[optionTable]
    okMask:(optionTable`status)=`OK;
    asc distinct optionTable[`snapshotDate] where okMask
 };

/ Multi-day replay: consecutive date pairs between startDate and endDate
.parser.barchart.multiDayReplay:{[optionTable;startDate;endDate;qty]
    allDates:.parser.barchart.availableDates optionTable;
    rangeDates:allDates where (allDates>=startDate) and allDates<=endDate;
    if[2>count rangeDates; '"fewer than 2 available dates in range"];
    allResults:();
    pairIdx:0;
    nPairs:(count rangeDates)-1;
    while[pairIdx<nPairs;
        entryDt:rangeDates pairIdx;
        exitDt:rangeDates pairIdx+1;
        oneReplay:.parser.barchart.__safeOneDayReplay[optionTable;entryDt;exitDt;qty];
        if[0<count oneReplay; allResults:allResults,oneReplay];
        pairIdx+:1];
    if[0=count allResults; '"no replay rows produced across date range"];
    allResults
 };

/ Aggregate marketPnl by exitDate
.parser.barchart.dailyPnl:{[replayTable]
    if[0=count replayTable; '"empty replay table"];
    exitDates:asc distinct replayTable`exitDate;
    resultRows:();
    cumPnl:0f;
    dIdx:0;
    while[dIdx<count exitDates;
        exitDt:exitDates dIdx;
        dayMask:(replayTable`exitDate)=exitDt;
        dayRows:replayTable where dayMask;
        statusCol:dayRows`status;
        okCnt:sum statusCol=`OK;
        errCnt:(count dayRows)-okCnt;
        okRows:dayRows where statusCol=`OK;
        dayPnl:$[0<count okRows;sum okRows`marketPnl;0f];
        cumPnl+:dayPnl;
        resultRows:resultRows,enlist `snapshotDate`dailyMarketPnl`cumulativeMarketPnl`tradeCount`okRows`errorRows`status`errorMessage!(
            exitDt;dayPnl;cumPnl;count dayRows;okCnt;errCnt;`OK;"");
        dIdx+:1];
    resultRows
 };

/ Aggregate marketPnl by contractId and exitDate
.parser.barchart.contractPnlSeries:{[replayTable]
    if[0=count replayTable; '"empty replay table"];
    contractIds:distinct replayTable`contractId;
    resultRows:();
    cIdx:0;
    while[cIdx<count contractIds;
        cid:contractIds cIdx;
        cidMask:(replayTable`contractId)=cid;
        cidRows:replayTable where cidMask;
        exitDates:asc distinct cidRows`exitDate;
        cumPnl:0f;
        eIdx:0;
        while[eIdx<count exitDates;
            exitDt:exitDates eIdx;
            dayMask:(cidRows`exitDate)=exitDt;
            dayRows:cidRows where dayMask;
            dayPnl:$[0<count dayRows;sum dayRows`marketPnl;0f];
            cumPnl+:dayPnl;
            resultRows:resultRows,enlist `contractId`exitDate`dailyMarketPnl`cumulativeMarketPnl`status`errorMessage!(
                cid;exitDt;dayPnl;cumPnl;`OK;"");
            eIdx+:1];
        cIdx+:1];
    resultRows
 };

/ Drawdown from daily PnL: running peak and loss from peak
.parser.barchart.drawdown:{[dailyPnlTable]
    if[0=count dailyPnlTable; '"empty dailyPnl table"];
    cumPnlCol:dailyPnlTable`cumulativeMarketPnl;
    snapshotCol:dailyPnlTable`snapshotDate;
    nRows:count dailyPnlTable;
    peakVals:nRows#0f;
    ddVals:nRows#0f;
    runPeak:0f;
    ddIdx:0;
    while[ddIdx<nRows;
        currentCum:cumPnlCol ddIdx;
        if[currentCum>runPeak; runPeak:currentCum];
        peakVals[ddIdx]:runPeak;
        ddVals[ddIdx]:runPeak-currentCum;
        ddIdx+:1];
    flip `snapshotDate`cumulativeMarketPnl`runningPeak`drawdown!(
        snapshotCol;cumPnlCol;peakVals;ddVals)
 };

/ Performance summary from daily PnL
.parser.barchart.performanceSummary:{[dailyPnlTable]
    if[0=count dailyPnlTable; '"empty dailyPnl table"];
    dailyPnlCol:dailyPnlTable`dailyMarketPnl;
    nObs:count dailyPnlCol;
    totalPnl:sum dailyPnlCol;
    meanPnl:avg dailyPnlCol;
    volPnl:dev dailyPnlCol;
    winDays:sum dailyPnlCol>0f;
    winRateVal:winDays%nObs;
    worstDay:min dailyPnlCol;
    bestDay:max dailyPnlCol;
    ddTable:.parser.barchart.drawdown dailyPnlTable;
    maxDd:max ddTable`drawdown;
    `observationCount`totalMarketPnl`meanDailyPnl`dailyPnlVolatility`winRate`worstDailyPnl`bestDailyPnl`maxDrawdown`status`errorMessage!(
        nObs;totalPnl;meanPnl;volPnl;winRateVal;worstDay;bestDay;maxDd;`OK;"")
 };

/ Overall backtest summary from replay table
.parser.barchart.backtestSummary:{[replayTable]
    if[0=count replayTable; '"empty replay table"];
    statusCol:replayTable`status;
    okCnt:sum statusCol=`OK;
    errCnt:(count replayTable)-okCnt;
    pnlCol:replayTable`marketPnl;
    contractIds:distinct replayTable`contractId;
    `rowCount`contractCount`startDate`endDate`totalMarketPnl`worstTradePnl`bestTradePnl`errorRows`status`errorMessage!(
        count replayTable;count contractIds;min replayTable`entryDate;max replayTable`exitDate;
        sum pnlCol;min pnlCol;max pnlCol;errCnt;`OK;"")
 };

/ ══════════════════════════════════════════════════════════════════
/ PRIVATE HELPERS
/ ══════════════════════════════════════════════════════════════════

/ Safe wrapper for oneDayReplay — returns empty list on error, logs warning
.parser.barchart.__safeOneDayReplay:{[optionTable;entryDt;exitDt;qty]
    @[{.parser.barchart.oneDayReplay[x 0;x 1;x 2;x 3]};(optionTable;entryDt;exitDt;qty);
        {[ed;xd;e] -1 "  replay skip ",string[ed],"->",string[xd],": ",e; ()}[entryDt;exitDt;]]
 };

-1 "backtest.q loaded — .parser.barchart backtest functions ready";
