/ commodityBacktest.q - commodity futures backtest (v0.32)

.backtest.commodity.longFrontRoll:{[futuresCurve;startDate;endDate;rollDaysBefore;qty]
    okMask:(futuresCurve`status)=`OK;
    dateMask:((futuresCurve`snapshotDate)>=startDate) and (futuresCurve`snapshotDate)<=endDate;
    curveRows:futuresCurve where okMask and dateMask;
    if[0=count curveRows; '"No curve data in date range"];
    availDates:asc distinct curveRows`snapshotDate;
    if[2>count availDates; '"Fewer than 2 dates in range"];
    resultRows:();
    currentTenor:1;
    prevPrice:0Nf;
    prevContract:`;
    cumPnl:0f;
    adIdx:0;
    while[adIdx<count availDates;
          sd:availDates adIdx;
          dayMask:(curveRows`snapshotDate)=sd;
          dayRows:curveRows where dayMask;
          / Get current tenor contract
          tenorMask:(dayRows`tenorRank)=currentTenor;
          tenorRows:dayRows where tenorMask;
          if[0=count tenorRows;
             / Tenor not available, try front
             tenorMask2:(dayRows`tenorRank)=1;
             tenorRows:dayRows where tenorMask2;
             currentTenor:1];
          if[0<count tenorRows;
             currentRow:tenorRows 0;
             currentPx:currentRow`settlementPrice;
             currentContract:currentRow`contractCode;
             mltVal:currentRow`contractMultiplier;
             dte:currentRow`daysToExpiry;
             / Check roll
             needsRoll:(dte<=rollDaysBefore) and currentTenor=1;
             / Get next contract for roll
             nextTenorMask:(dayRows`tenorRank)=currentTenor+1;
             nextRows:dayRows where nextTenorMask;
             nextContract:$[0<count nextRows;(nextRows 0)`contractCode;`];
             / Daily PnL
             dayPnl:$[(not null prevPrice) and prevContract=currentContract;
                      qty*mltVal*currentPx-prevPrice;0f];
             cumPnl+:dayPnl;
             rollPnl:0f;
             rolledFlag:0b;
             actionSym:`hold;
             if[adIdx=0; actionSym:`enter];
             if[needsRoll and 0<count nextRows;
                actionSym:`roll;
                rolledFlag:1b;
                rollPnl:dayPnl;
                / Switch to next contract: update prev so next day tracks new contract
                prevPrice:(nextRows 0)`settlementPrice;
                prevContract:(nextRows 0)`contractCode;
                currentTenor:currentTenor+1];
             resultRows:resultRows,enlist `snapshotDate`commodity`currentContract`nextContract`action`quantity`contractMultiplier`entryPrice`exitPrice`dailyPnl`cumulativePnl`rolled`rollPnl`status`errorMessage!(
                 sd;currentRow`commodity;currentContract;nextContract;actionSym;qty;mltVal;
                 prevPrice;currentPx;dayPnl;cumPnl;rolledFlag;rollPnl;`OK;"");
             if[not rolledFlag;
                prevPrice:currentPx;
                prevContract:currentContract]];
          adIdx+:1];
    resultRows
 };

.backtest.commodity.calendarSpreadReplay:{[futuresCurve;nearTenor;farTenor;startDate;endDate;qty]
    okMask:(futuresCurve`status)=`OK;
    dateMask:((futuresCurve`snapshotDate)>=startDate) and (futuresCurve`snapshotDate)<=endDate;
    curveRows:futuresCurve where okMask and dateMask;
    availDates:asc distinct curveRows`snapshotDate;
    if[2>count availDates; '"Fewer than 2 dates for spread replay"];
    resultRows:();
    prevSpread:0Nf;
    cumPnl:0f;
    csIdx:0;
    while[csIdx<count availDates;
          sd:availDates csIdx;
          dayMask:(curveRows`snapshotDate)=sd;
          dayRows:curveRows where dayMask;
          nearMask:(dayRows`tenorRank)=nearTenor;
          farMask:(dayRows`tenorRank)=farTenor;
          nearRows:dayRows where nearMask;
          farRows:dayRows where farMask;
          if[(0<count nearRows) and 0<count farRows;
             nearPx:(nearRows 0)`settlementPrice;
             farPx:(farRows 0)`settlementPrice;
             spreadVal:nearPx-farPx;
             mltVal:(nearRows 0)`contractMultiplier;
             dayPnl:$[not null prevSpread;qty*mltVal*spreadVal-prevSpread;0f];
             cumPnl+:dayPnl;
             resultRows:resultRows,enlist `snapshotDate`nearContract`farContract`nearPrice`farPrice`spread`dailyPnl`cumulativePnl`status`errorMessage!(
                 sd;(nearRows 0)`contractCode;(farRows 0)`contractCode;nearPx;farPx;spreadVal;dayPnl;cumPnl;`OK;"");
             prevSpread:spreadVal];
          csIdx+:1];
    resultRows
 };

.backtest.commodity.dailyPnl:{[replayTable]
    if[0=count replayTable; '"Empty replay table"];
    exitDates:asc distinct replayTable`snapshotDate;
    resultRows:();
    cumPnl:0f;
    dpIdx:0;
    while[dpIdx<count exitDates;
          sd:exitDates dpIdx;
          dayMask:(replayTable`snapshotDate)=sd;
          dayRows:replayTable where dayMask;
          dayPnl:sum dayRows`dailyPnl;
          cumPnl+:dayPnl;
          resultRows:resultRows,enlist `snapshotDate`dailyPnl`cumulativePnl`status`errorMessage!(sd;dayPnl;cumPnl;`OK;"");
          dpIdx+:1];
    resultRows
 };

.backtest.commodity.performanceSummary:{[dailyPnlTable]
    if[0=count dailyPnlTable; '"Empty daily PnL table"];
    pnlCol:dailyPnlTable`dailyPnl;
    nObs:count pnlCol;
    cumCol:dailyPnlTable`cumulativePnl;
    / Max drawdown
    runPeak:0f;
    maxDd:0f;
    ddIdx:0;
    while[ddIdx<nObs;
          if[cumCol[ddIdx]>runPeak; runPeak:cumCol ddIdx];
          dd:runPeak-cumCol ddIdx;
          if[dd>maxDd; maxDd:dd];
          ddIdx+:1];
    `observationCount`totalPnl`meanDailyPnl`dailyPnlVolatility`winRate`worstDailyPnl`bestDailyPnl`maxDrawdown`status`errorMessage!(
        nObs;sum pnlCol;avg pnlCol;dev pnlCol;(sum pnlCol>0f)%nObs;min pnlCol;max pnlCol;maxDd;`OK;"")
 };
