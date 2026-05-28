/ regression.q - regression checks: PV changes, large PnL, unexplained PnL

/ --- Public ---

.regression.comparePricingRuns:{[currentPricing;previousPricing;toleranceVal]
    numTrades:count currentPricing;
    resultList:();
    loopIdx:0;
    while[loopIdx<numTrades;
        currentRow:currentPricing loopIdx;
        comparisonRow:.regression.__buildComparisonRow[currentRow;previousPricing;toleranceVal];
        resultList:resultList,enlist comparisonRow;
        loopIdx+:1];
    resultList
 };

.regression.flagLargePnL:{[pnlExplainResult;absoluteThreshold]
    okRows:pnlExplainResult where pnlExplainResult[`status]=`OK;
    okRows where (abs okRows`actualPnL)>absoluteThreshold
 };

.regression.flagLargeUnexplainedPnL:{[pnlExplainResult;absoluteThreshold]
    okRows:pnlExplainResult where pnlExplainResult[`status]=`OK;
    okRows where (abs okRows`unexplainedPnL)>absoluteThreshold
 };

/ --- Internal ---

.regression.__buildComparisonRow:{[currentRow;previousPricing;toleranceVal]
    tradeIdVal:currentRow`tradeId;
    underlyingVal:currentRow`underlying;
    currentPV:currentRow`unitPrice;
    prevRows:previousPricing where previousPricing[`tradeId]=tradeIdVal;
    if[0=count prevRows;
        :`tradeId`underlyingSym`previousPV`currentPV`pvChange`pvChangePct`status`warningMessage!(
            tradeIdVal;underlyingVal;0Nf;currentPV;0Nf;0Nf;`NO_PREVIOUS;"No previous pricing found")];
    previousPV:(prevRows`unitPrice)[0];
    pvChange:currentPV-previousPV;
    pvChangePct:0f;
    if[not previousPV=0f; pvChangePct:pvChange%previousPV];
    statusVal:`OK;
    warningMsg:"";
    if[(abs pvChange)>toleranceVal;
        statusVal:`WARNING;
        warningMsg:"PV change ",string[pvChange]," exceeds tolerance ",string toleranceVal];
    `tradeId`underlyingSym`previousPV`currentPV`pvChange`pvChangePct`status`warningMessage!(
        tradeIdVal;underlyingVal;previousPV;currentPV;pvChange;pvChangePct;statusVal;warningMsg)
 };
