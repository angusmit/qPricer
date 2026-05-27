/ dashboard.q - risk dashboard summary (v0.30)

.dashboard.__safeGet:{[dict;keyName;defaultVal]
    if[not 99h=type dict; :defaultVal];
    if[0=count dict; :defaultVal];
    if[keyName in key dict; :dict keyName];
    defaultVal
 };

.dashboard.__safeSum:{[tbl;colName]
    if[0=count tbl; :0Nf];
    if[not colName in key tbl 0; :0Nf];
    sum tbl colName
 };

.dashboard.portfolioTile:{[pricingResult]
    if[0=count pricingResult; :`tradeCount`pricedOkCount`pricingErrorCount`totalPV!(0;0;0;0Nf)];
    totalTrades:count pricingResult;
    okMask:$[`status in key pricingResult 0;(pricingResult`status)=`OK;totalTrades#1b];
    okCnt:sum okMask;
    totalPV:.dashboard.__safeSum[pricingResult;`notionalPrice];
    `tradeCount`pricedOkCount`pricingErrorCount`totalPV!(totalTrades;okCnt;totalTrades-okCnt;totalPV)
 };

.dashboard.greeksTile:{[greekResult]
    if[0=count greekResult; :`totalDelta`totalGamma`totalVega`totalTheta`totalRho!(0Nf;0Nf;0Nf;0Nf;0Nf)];
    deltaVal:.dashboard.__safeSum[greekResult;`delta];
    gammaVal:.dashboard.__safeSum[greekResult;`gamma];
    vegaVal:.dashboard.__safeSum[greekResult;`vega];
    thetaVal:.dashboard.__safeSum[greekResult;`theta];
    rhoVal:.dashboard.__safeSum[greekResult;`rho];
    `totalDelta`totalGamma`totalVega`totalTheta`totalRho!(deltaVal;gammaVal;vegaVal;thetaVal;rhoVal)
 };

.dashboard.varTile:{[varReport]
    if[0=count varReport; :`var95`es95!(0Nf;0Nf)];
    firstRow:varReport 0;
    `var95`es95!(firstRow`valueAtRisk;firstRow`expectedShortfall)
 };

.dashboard.historicalTile:{[histReplayResult;histVarReport;worstEvents]
    histVar95:$[0<count histVarReport;(histVarReport 0)`valueAtRisk;0Nf];
    histEs95:$[0<count histVarReport;(histVarReport 0)`expectedShortfall;0Nf];
    worstLoss:$[0<count worstEvents;(worstEvents 0)`loss;0Nf];
    `historicalVar95`historicalEs95`worstHistoricalLoss!(histVar95;histEs95;worstLoss)
 };

.dashboard.limitTile:{[limitDash]
    `limitOkCount`limitWarningCount`limitBreachCount`limitErrorCount!(
        .dashboard.__safeGet[limitDash;`okCount;0];
        .dashboard.__safeGet[limitDash;`warningCount;0];
        .dashboard.__safeGet[limitDash;`breachCount;0];
        .dashboard.__safeGet[limitDash;`errorCount;0])
 };

.dashboard.modelCheckTile:{[mcSummary]
    `modelCheckPassedCount`modelCheckFailedCount!(
        .dashboard.__safeGet[mcSummary;`passedCount;0];
        .dashboard.__safeGet[mcSummary;`failedCount;0])
 };

.dashboard.dashboardSummary:{[dailyRiskResult]
    pTile:.dashboard.portfolioTile[dailyRiskResult`pricingResult];
    gTile:.dashboard.greeksTile[dailyRiskResult`greekResult];
    vTile:.dashboard.varTile[dailyRiskResult`varReport];
    hTile:.dashboard.historicalTile[dailyRiskResult`historicalReplayResult;dailyRiskResult`historicalVarReport;dailyRiskResult`worstHistoricalEvents];
    lTile:.dashboard.limitTile[dailyRiskResult`limitDashboard];
    mTile:.dashboard.modelCheckTile[dailyRiskResult`modelCheckSummary];
    / Overall status
    overallStat:`OK;
    if[lTile[`limitBreachCount]>0; overallStat:`breach];
    if[(overallStat=`OK) and lTile[`limitWarningCount]>0; overallStat:`warning];
    if[pTile[`pricingErrorCount]>0; overallStat:`error];
    if[mTile[`modelCheckFailedCount]>0; if[overallStat=`OK; overallStat:`warning]];
    / Build result dict step by step to avoid right-to-left chain type issues
    result:`runLabel`valuationDate!("DailyRisk";.z.D);
    result[`tradeCount]:pTile`tradeCount;
    result[`pricedOkCount]:pTile`pricedOkCount;
    result[`pricingErrorCount]:pTile`pricingErrorCount;
    result[`totalPV]:pTile`totalPV;
    result[`totalDelta]:gTile`totalDelta;
    result[`totalGamma]:gTile`totalGamma;
    result[`totalVega]:gTile`totalVega;
    result[`totalTheta]:gTile`totalTheta;
    result[`totalRho]:gTile`totalRho;
    result[`var95]:vTile`var95;
    result[`es95]:vTile`es95;
    result[`historicalVar95]:hTile`historicalVar95;
    result[`historicalEs95]:hTile`historicalEs95;
    result[`worstHistoricalLoss]:hTile`worstHistoricalLoss;
    result[`limitOkCount]:lTile`limitOkCount;
    result[`limitWarningCount]:lTile`limitWarningCount;
    result[`limitBreachCount]:lTile`limitBreachCount;
    result[`limitErrorCount]:lTile`limitErrorCount;
    result[`modelCheckPassedCount]:mTile`modelCheckPassedCount;
    result[`modelCheckFailedCount]:mTile`modelCheckFailedCount;
    result[`overallStatus]:overallStat;
    result[`status]:`OK;
    result[`errorMessage]:"";
    result
 };
