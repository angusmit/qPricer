/ backtest.q - Barchart multi-day replay & performance analytics
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

/ Safe wrapper for oneDayReplay - returns empty list on error, logs warning
.parser.barchart.__safeOneDayReplay:{[optionTable;entryDt;exitDt;qty]
    @[{.parser.barchart.oneDayReplay[x 0;x 1;x 2;x 3]};(optionTable;entryDt;exitDt;qty);
      {[ed;xd;e] -1 "  replay skip ",string[ed],"->",string[xd],": ",e; ()}[entryDt;exitDt;]]
 };

-1 "backtest.q loaded - .parser.barchart backtest functions ready";

/ ══════════════════════════════════════════════════════════════════
/ v0.31.1 - MODEL VS MARKET COMPARISON
/ ══════════════════════════════════════════════════════════════════

/ ── Private: analytical BS pricing & Greeks ───────────────────

.parser.barchart.__normalPdf:{[x] exp[neg 0.5*x*x]%sqrt 2f*acos neg 1f};

/ Price one option row using existing qFDM BS closed-form
.parser.barchart.__bsPrice:{[optType;spotVal;strikeVal;tauVal;rateVal;divYVal;volVal]
    @[{.validation.blackScholesClosedForm[x 0;x 1;x 2;x 3;x 4;x 5;x 6]};(optType;spotVal;strikeVal;tauVal;rateVal;divYVal;volVal);{0Nf}]
 };

/ Analytical BS Greeks for one row
.parser.barchart.__bsGreeks:{[optType;spotVal;strikeVal;tauVal;rateVal;divYVal;volVal]
    sqrtT:sqrt tauVal;
    d1Val:.validation.__d1[spotVal;strikeVal;tauVal;rateVal;divYVal;volVal];
    d2Val:.validation.__d2[d1Val;volVal;tauVal];
    nd1:.validation.__normalCdf d1Val;
    nd2:.validation.__normalCdf d2Val;
    nNd1:.validation.__normalCdf neg d1Val;
    nNd2:.validation.__normalCdf neg d2Val;
    pdfD1:.parser.barchart.__normalPdf d1Val;
    eqT:exp neg divYVal*tauVal;
    erT:exp neg rateVal*tauVal;
    / Delta
    deltaVal:$[optType=`call;eqT*nd1;neg eqT*nNd1];
    / Gamma
    gammaVal:eqT*pdfD1%(spotVal*volVal*sqrtT);
    / Vega (per 1% = divide by 100)
    vegaVal:spotVal*eqT*pdfD1*sqrtT%100f;
    / Theta (per day = divide by 365)
    thetaPart1:neg spotVal*eqT*pdfD1*volVal%(2f*sqrtT);
    thetaVal:$[optType=`call;
               (thetaPart1-(rateVal*strikeVal*erT*nd2))+(divYVal*spotVal*eqT*nd1);
               (thetaPart1+(rateVal*strikeVal*erT*nNd2))-(divYVal*spotVal*eqT*nNd1)];
    thetaVal:thetaVal%365f;
    / Rho (per 1% = divide by 100)
    rhoVal:$[optType=`call;
             strikeVal*tauVal*erT*nd2%100f;
             neg strikeVal*tauVal*erT*nNd2%100f];
    `modelDelta`modelGamma`modelVega`modelTheta`modelRho!(deltaVal;gammaVal;vegaVal;thetaVal;rhoVal)
 };

/ ── Public: price all option rows ─────────────────────────────

.parser.barchart.priceOptionRows:{[optionTable;pricingConfig]
    rateVal:pricingConfig`riskFreeRate;
    divYVal:pricingConfig`dividendYield;
    nRows:count optionTable;
    / Pre-allocate result columns
    modelPrices:nRows#0Nf;
    modelErrors:nRows#0Nf;
    vendorTheoErrors:nRows#0Nf;
    mDeltas:nRows#0Nf;  mGammas:nRows#0Nf;  mVegas:nRows#0Nf;  mThetas:nRows#0Nf;  mRhos:nRows#0Nf;
    deltaErrs:nRows#0Nf; gammaErrs:nRows#0Nf; vegaErrs:nRows#0Nf; thetaErrs:nRows#0Nf; rhoErrs:nRows#0Nf;
    pStatuses:nRows#`OK;
    pErrors:nRows#enlist "";
    rIdx:0;
    while[rIdx<nRows;
          rowData:optionTable rIdx;
          canPrice:rowData[`status]=`OK;
          if[not canPrice;
             pStatuses[rIdx]:`skipped;
             pErrors[rIdx]:"source row not OK"];
          if[canPrice;
             tauVal:rowData`tau;
             volVal:rowData`impliedVolatility;
             spotVal:rowData`spot;
             optType:rowData`optionType;
             strikeVal:rowData`strike;
             if[(null tauVal) or tauVal<=0f; pStatuses[rIdx]:`error; pErrors[rIdx]:"non-positive tau"; canPrice:0b];
             if[canPrice and ((null volVal) or volVal<=0f); pStatuses[rIdx]:`error; pErrors[rIdx]:"non-positive vol"; canPrice:0b];
             if[canPrice and ((null spotVal) or spotVal<=0f); pStatuses[rIdx]:`error; pErrors[rIdx]:"non-positive spot"; canPrice:0b]];
          if[canPrice;
             mp:.parser.barchart.__bsPrice[optType;spotVal;strikeVal;tauVal;rateVal;divYVal;volVal];
             if[null mp; pStatuses[rIdx]:`error; pErrors[rIdx]:"BS pricing failed"; canPrice:0b]];
          if[canPrice;
             modelPrices[rIdx]:mp;
             modelErrors[rIdx]:mp-rowData`marketPrice;
             if[not null rowData`vendorTheo; vendorTheoErrors[rIdx]:mp-rowData`vendorTheo];
             greeks:@[{.parser.barchart.__bsGreeks[x 0;x 1;x 2;x 3;x 4;x 5;x 6]};(optType;spotVal;strikeVal;tauVal;rateVal;divYVal;volVal);{`modelDelta`modelGamma`modelVega`modelTheta`modelRho!(0Nf;0Nf;0Nf;0Nf;0Nf)}];
             mDeltas[rIdx]:greeks`modelDelta;
             mGammas[rIdx]:greeks`modelGamma;
             mVegas[rIdx]:greeks`modelVega;
             mThetas[rIdx]:greeks`modelTheta;
             mRhos[rIdx]:greeks`modelRho;
             deltaErrs[rIdx]:greeks[`modelDelta]-rowData`vendorDelta;
             gammaErrs[rIdx]:greeks[`modelGamma]-rowData`vendorGamma;
             vegaErrs[rIdx]:greeks[`modelVega]-rowData`vendorVega;
             thetaErrs[rIdx]:greeks[`modelTheta]-rowData`vendorTheta;
             rhoErrs[rIdx]:greeks[`modelRho]-rowData`vendorRho];
          rIdx+:1];
    / Append new columns to optionTable
    optionTable,'flip `modelPrice`modelError`vendorTheoError`modelDelta`modelGamma`modelVega`modelTheta`modelRho`deltaError`gammaError`vegaError`thetaError`rhoError`pricingStatus`pricingErrorMessage!(
        modelPrices;modelErrors;vendorTheoErrors;mDeltas;mGammas;mVegas;mThetas;mRhos;
        deltaErrs;gammaErrs;vegaErrs;thetaErrs;rhoErrs;pStatuses;pErrors)
 };

/ ── Public: model replay (join model prices with market replay) ─

.parser.barchart.modelReplay:{[replayTable;pricedOptionTable]
    / Build lookup: (contractId, snapshotDate) -> row index
    cidCol:pricedOptionTable`contractId;
    dateCol:pricedOptionTable`snapshotDate;
    lookupKeys:{`$string[x],"|",string y}'[cidCol;dateCol];
    lookupDict:lookupKeys!til count pricedOptionTable;
    nRows:count replayTable;
    entryModelPrices:nRows#0Nf;  exitModelPrices:nRows#0Nf;
    modelPnls:nRows#0Nf;         pnlDiffs:nRows#0Nf;
    entryModelErrs:nRows#0Nf;    exitModelErrs:nRows#0Nf;
    entryVTheos:nRows#0Nf;       exitVTheos:nRows#0Nf;
    entryVTheoErrs:nRows#0Nf;    exitVTheoErrs:nRows#0Nf;
    pStatuses:nRows#`OK;
    pErrors:nRows#enlist "";
    rIdx:0;
    while[rIdx<nRows;
          rRow:replayTable rIdx;
          cid:rRow`contractId;
          entryKey:`$string[cid],"|",string rRow`entryDate;
          exitKey:`$string[cid],"|",string rRow`exitDate;
          hasEntry:entryKey in lookupKeys;
          hasExit:exitKey in lookupKeys;
          if[hasEntry and hasExit;
             eIdx:lookupDict entryKey;
             xIdx:lookupDict exitKey;
             eRow:pricedOptionTable eIdx;
             xRow:pricedOptionTable xIdx;
             emp:eRow`modelPrice;
             xmp:xRow`modelPrice;
             entryModelPrices[rIdx]:emp;
             exitModelPrices[rIdx]:xmp;
             if[(not null emp) and not null xmp;
                mPnl:rRow[`quantity]*.parser.barchart.cfg.contractMult*xmp-emp;
                modelPnls[rIdx]:mPnl;
                pnlDiffs[rIdx]:mPnl-rRow`marketPnl];
             entryModelErrs[rIdx]:emp-rRow`entryMarketPrice;
             exitModelErrs[rIdx]:xmp-rRow`exitMarketPrice;
             if[`vendorTheo in key eRow;
                entryVTheos[rIdx]:eRow`vendorTheo;
                exitVTheos[rIdx]:xRow`vendorTheo;
                if[not null eRow`vendorTheo; entryVTheoErrs[rIdx]:emp-eRow`vendorTheo];
                if[not null xRow`vendorTheo; exitVTheoErrs[rIdx]:xmp-xRow`vendorTheo]]];
          if[not hasEntry and hasExit;
             pStatuses[rIdx]:`error; pErrors[rIdx]:"missing entry model price"];
          if[hasEntry and not hasExit;
             pStatuses[rIdx]:`error; pErrors[rIdx]:"missing exit model price"];
          if[(not hasEntry) and not hasExit;
             pStatuses[rIdx]:`error; pErrors[rIdx]:"missing entry and exit model prices"];
          rIdx+:1];
    replayTable,'flip `entryModelPrice`exitModelPrice`modelPnl`pnlDifference`entryModelError`exitModelError`entryVendorTheo`exitVendorTheo`entryVendorTheoError`exitVendorTheoError`pricingStatus`pricingErrorMessage!(
        entryModelPrices;exitModelPrices;modelPnls;pnlDiffs;entryModelErrs;exitModelErrs;
        entryVTheos;exitVTheos;entryVTheoErrs;exitVTheoErrs;pStatuses;pErrors)
 };

/ ── Public: model comparison summary ──────────────────────────

.parser.barchart.modelComparisonSummary:{[modelReplayTable]
    if[0=count modelReplayTable; '"empty model replay table"];
    statusCol:modelReplayTable`pricingStatus;
    okMask:statusCol=`OK;
    okRows:modelReplayTable where okMask;
    okCnt:count okRows;
    errCnt:(count modelReplayTable)-okCnt;
    mktPnls:okRows`marketPnl;
    mdlPnls:okRows`modelPnl;
    pnlDiffs:okRows`pnlDifference;
    entryErrs:okRows`entryModelError;
    exitErrs:okRows`exitModelError;
    `rowCount`okRows`errorRows`totalMarketPnl`totalModelPnl`totalPnlDifference`meanEntryModelError`meanExitModelError`meanAbsEntryModelError`meanAbsExitModelError`meanAbsPnlDifference`worstPnlDifference`bestPnlDifference`status`errorMessage!(
        count modelReplayTable;okCnt;errCnt;
        sum mktPnls;sum mdlPnls;sum pnlDiffs;
        avg entryErrs;avg exitErrs;avg abs entryErrs;avg abs exitErrs;avg abs pnlDiffs;
        min pnlDiffs;max pnlDiffs;`OK;"")
 };

/ ── Public: error bucket analysis ─────────────────────────────

/ Helper for bucket aggregation
.parser.barchart.__bucketStats:{[bucketName;bucketRows]
    modelErrs:bucketRows`modelError;
    validME:modelErrs where not null modelErrs;
    vtErrs:bucketRows`vendorTheoError;
    validVT:vtErrs where not null vtErrs;
    `bucket`rowCount`meanModelError`meanAbsModelError`meanVendorTheoError`meanAbsVendorTheoError!(
        bucketName;count bucketRows;
        $[0<count validME;avg validME;0Nf];
        $[0<count validME;avg abs validME;0Nf];
        $[0<count validVT;avg validVT;0Nf];
        $[0<count validVT;avg abs validVT;0Nf])
 };

.parser.barchart.errorByMoneyness:{[pricedOptionTable]
    okMask:(pricedOptionTable`pricingStatus)=`OK;
    okRows:pricedOptionTable where okMask;
    moneyCol:okRows`moneyness;
    lowMask:moneyCol<0.9;
    atmMask:(moneyCol>=0.9) and moneyCol<=1.1;
    highMask:moneyCol>1.1;
    resultRows:();
    resultRows:resultRows,enlist .parser.barchart.__bucketStats[`lowStrike;okRows where lowMask];
    resultRows:resultRows,enlist .parser.barchart.__bucketStats[`atm;okRows where atmMask];
    resultRows:resultRows,enlist .parser.barchart.__bucketStats[`highStrike;okRows where highMask];
    resultRows
 };

.parser.barchart.errorByDte:{[pricedOptionTable]
    okMask:(pricedOptionTable`pricingStatus)=`OK;
    okRows:pricedOptionTable where okMask;
    dteCol:okRows`dte;
    shortMask:dteCol<=30i;
    medMask:(dteCol>30i) and dteCol<=60i;
    longMask:dteCol>60i;
    resultRows:();
    resultRows:resultRows,enlist .parser.barchart.__bucketStats[`short;okRows where shortMask];
    resultRows:resultRows,enlist .parser.barchart.__bucketStats[`medium;okRows where medMask];
    resultRows:resultRows,enlist .parser.barchart.__bucketStats[`long;okRows where longMask];
    resultRows
 };

.parser.barchart.errorByOptionType:{[pricedOptionTable]
    okMask:(pricedOptionTable`pricingStatus)=`OK;
    okRows:pricedOptionTable where okMask;
    optTypes:distinct okRows`optionType;
    resultRows:();
    oIdx:0;
    while[oIdx<count optTypes;
          ot:optTypes oIdx;
          otRows:okRows where (okRows`optionType)=ot;
          resultRows:resultRows,enlist .parser.barchart.__bucketStats[ot;otRows];
          oIdx+:1];
    resultRows
 };
