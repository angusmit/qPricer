/ limits.q - risk limit monitoring (v0.29)

.limits.__requiredLimitCols:`limitId`scopeType`scopeValue`metricName`limitValue`warningPct`hardLimitPct`direction`enabled;
.limits.__supportedDirections:`lessThan`greaterThan`absLessThan;

.limits.validateLimitTable:{[limitTable]
    if[0=count limitTable; '"Limit table is empty"];
    firstRow:limitTable 0;
    tableCols:key firstRow;
    missingCols:.limits.__requiredLimitCols where not .limits.__requiredLimitCols in tableCols;
    if[0<count missingCols; '"Missing limit columns: ",raze ", " ,/: string missingCols];
    limitIds:limitTable`limitId;
    if[(count limitIds)<>count distinct limitIds; '"Duplicate limitId values"];
    limitVals:limitTable`limitValue;
    if[any limitVals<=0f; '"limitValue must be positive"];
    warnPcts:limitTable`warningPct;
    hardPcts:limitTable`hardLimitPct;
    if[any hardPcts<=0f; '"hardLimitPct must be positive"];
    if[any warnPcts>hardPcts; '"warningPct must be <= hardLimitPct"];
    directions:limitTable`direction;
    badDirs:directions where not directions in .limits.__supportedDirections;
    if[0<count badDirs; '"Invalid direction: ",string first badDirs];
 };

.limits.metricValue:{[metricTable;metricName;scopeType;scopeValue]
    nameMask:(metricTable`metricName)=metricName;
    scopeMask:(metricTable`scopeType)=scopeType;
    valMask:(metricTable`scopeValue)=scopeValue;
    allMask:nameMask and scopeMask and valMask;
    matchRows:metricTable where allMask;
    if[0=count matchRows; :0Nf];
    (matchRows 0)`metricValue
 };

.limits.evaluateLimit:{[limitRow;metricVal]
    limitId:limitRow`limitId;
    scopeType:limitRow`scopeType;
    scopeValue:limitRow`scopeValue;
    metricName:limitRow`metricName;
    limitVal:limitRow`limitValue;
    warnPct:limitRow`warningPct;
    hardPct:limitRow`hardLimitPct;
    directionVal:limitRow`direction;
    warnThreshold:warnPct*limitVal;
    hardThreshold:hardPct*limitVal;
    if[null metricVal;
        :`limitId`scopeType`scopeValue`metricName`metricValue`limitValue`warningThreshold`hardThreshold`utilisation`severity`passed`status`errorMessage!(
            limitId;scopeType;scopeValue;metricName;0Nf;limitVal;warnThreshold;hardThreshold;0Nf;`error;0b;`ERROR;"Metric value not found")];
    / Compute monitored value based on direction
    monitoredVal:$[directionVal=`absLessThan;abs metricVal;metricVal];
    / Compute utilisation
    utilisationVal:$[hardThreshold>0f;monitoredVal%hardThreshold;0Nf];
    / Evaluate severity
    severityVal:`OK;
    passedVal:1b;
    $[directionVal in `lessThan`absLessThan;[
        if[monitoredVal>hardThreshold; severityVal:`breach; passedVal:0b];
        if[(severityVal=`OK) and monitoredVal>warnThreshold; severityVal:`warning; passedVal:1b]];
      directionVal=`greaterThan;[
        if[monitoredVal<hardThreshold; severityVal:`breach; passedVal:0b];
        if[(severityVal=`OK) and monitoredVal<warnThreshold; severityVal:`warning; passedVal:1b];
        utilisationVal:$[monitoredVal>0f;hardThreshold%monitoredVal;0Nf]]];
    `limitId`scopeType`scopeValue`metricName`metricValue`limitValue`warningThreshold`hardThreshold`utilisation`severity`passed`status`errorMessage!(
        limitId;scopeType;scopeValue;metricName;metricVal;limitVal;warnThreshold;hardThreshold;utilisationVal;severityVal;passedVal;`OK;"")
 };

.limits.evaluateLimits:{[limitTable;metricTable]
    .limits.validateLimitTable limitTable;
    if[0=count metricTable; '"Metric table is empty"];
    resultRows:();
    lIdx:0;
    while[lIdx<count limitTable;
        limitRow:limitTable lIdx;
        if[limitRow`enabled;
            metricVal:.limits.metricValue[metricTable;limitRow`metricName;limitRow`scopeType;limitRow`scopeValue];
            resultRows:resultRows,enlist .limits.evaluateLimit[limitRow;metricVal]];
        lIdx+:1];
    resultRows
 };

.limits.buildMetricTable:{[varReport;greekReport;pnlReport;historicalReport]
    resultRows:();
    / VaR report
    if[not (::)~varReport;
        vIdx:0;
        while[vIdx<count varReport;
            vRow:varReport vIdx;
            clName:`$"VaR",string `int$100*vRow`confidenceLevel;
            resultRows:resultRows,enlist `scopeType`scopeValue`metricName`metricValue`source!(`portfolio;`ALL;clName;vRow`valueAtRisk;`var);
            esName:`$"ES",string `int$100*vRow`confidenceLevel;
            resultRows:resultRows,enlist `scopeType`scopeValue`metricName`metricValue`source!(`portfolio;`ALL;esName;vRow`expectedShortfall;`var);
            vIdx+:1]];
    / Greek report
    if[not (::)~greekReport;
        gIdx:0;
        while[gIdx<count greekReport;
            gRow:greekReport gIdx;
            scopeType:$[`scopeType in key gRow;gRow`scopeType;`portfolio];
            scopeValue:$[`scopeValue in key gRow;gRow`scopeValue;`ALL];
            greekCols:`DeltaCash`GammaCash`VegaCash`ThetaCash`RhoCash;
            gcIdx:0;
            while[gcIdx<count greekCols;
                gName:greekCols gcIdx;
                if[gName in key gRow;
                    resultRows:resultRows,enlist `scopeType`scopeValue`metricName`metricValue`source!(scopeType;scopeValue;gName;gRow gName;`greeks)];
                gcIdx+:1];
            gIdx+:1]];
    / PnL/worst-loss report
    if[not (::)~pnlReport;
        pIdx:0;
        while[pIdx<count pnlReport;
            pRow:pnlReport pIdx;
            if[`loss in key pRow;
                resultRows:resultRows,enlist `scopeType`scopeValue`metricName`metricValue`source!(`portfolio;`ALL;`WorstLoss;pRow`loss;`pnl)];
            pIdx+:1]];
    / Historical replay report
    if[not (::)~historicalReport;
        hIdx:0;
        while[hIdx<count historicalReport;
            hRow:historicalReport hIdx;
            if[`loss in key hRow;
                resultRows:resultRows,enlist `scopeType`scopeValue`metricName`metricValue`source!(`portfolio;`ALL;`HistWorstLoss;hRow`loss;`historical)];
            hIdx+:1]];
    resultRows
 };

.limits.checkVarLimits:{[varReport;limitTable]
    metricTable:.limits.buildMetricTable[varReport;(::);(::);(::)];
    .limits.evaluateLimits[limitTable;metricTable]
 };

.limits.checkGreekLimits:{[greekReport;limitTable]
    metricTable:.limits.buildMetricTable[(::);greekReport;(::);(::)];
    .limits.evaluateLimits[limitTable;metricTable]
 };

.limits.checkPnlLimits:{[pnlReport;limitTable]
    metricTable:.limits.buildMetricTable[(::);(::);pnlReport;(::)];
    .limits.evaluateLimits[limitTable;metricTable]
 };

.limits.checkHistoricalLimits:{[historicalReport;limitTable]
    metricTable:.limits.buildMetricTable[(::);(::);(::);historicalReport];
    .limits.evaluateLimits[limitTable;metricTable]
 };

.limits.limitBreaches:{[limitCheckResult]
    sevCol:limitCheckResult`severity;
    limitCheckResult where sevCol=`breach
 };

.limits.limitWarnings:{[limitCheckResult]
    sevCol:limitCheckResult`severity;
    limitCheckResult where sevCol=`warning
 };
