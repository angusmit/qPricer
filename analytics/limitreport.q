/ limitreport.q - risk limit reporting (v0.29)

.limitreport.breachReport:{[limitCheckResult]
    .limits.limitBreaches limitCheckResult
 };

.limitreport.warningReport:{[limitCheckResult]
    .limits.limitWarnings limitCheckResult
 };

.limitreport.summaryBySeverity:{[limitCheckResult]
    allSeverities:`OK`warning`breach`error;
    sevCol:limitCheckResult`severity;
    utilisationCol:limitCheckResult`utilisation;
    resultRows:();
    sIdx:0;
    while[sIdx<count allSeverities;
        sev:allSeverities sIdx;
        sevMask:sevCol=sev;
        sevCount:sum sevMask;
        sevUtils:utilisationCol where sevMask;
        maxUtil:$[0<count sevUtils;max sevUtils;0Nf];
        resultRows:resultRows,enlist `severity`limitCount`maxUtilisation!(sev;sevCount;maxUtil);
        sIdx+:1];
    resultRows
 };

.limitreport.summaryByScope:{[limitCheckResult]
    scopeTypes:limitCheckResult`scopeType;
    scopeValues:limitCheckResult`scopeValue;
    seen:();
    resultRows:();
    rIdx:0;
    while[rIdx<count limitCheckResult;
        st:scopeTypes rIdx;
        scopeVal:scopeValues rIdx;
        pairKey:(st;scopeVal);
        isNew:$[0=count seen;1b;not pairKey in seen];
        if[isNew;
            seen:seen,enlist pairKey;
            stMask:scopeTypes=st;
            svMask:scopeValues=scopeVal;
            bothMask:stMask and svMask;
            scopeRows:limitCheckResult where bothMask;
            sevCol:scopeRows`severity;
            breachCnt:sum sevCol=`breach;
            warnCnt:sum sevCol=`warning;
            utilCol:scopeRows`utilisation;
            validUtils:utilCol where not null utilCol;
            maxUtil:$[0<count validUtils;max validUtils;0Nf];
            worstSev:`OK;
            if[0<sum sevCol=`error; worstSev:`error];
            if[0<sum sevCol=`breach; worstSev:`breach];
            if[(worstSev=`OK) and 0<sum sevCol=`warning; worstSev:`warning];
            resultRows:resultRows,enlist `scopeType`scopeValue`limitCount`breachCount`warningCount`maxUtilisation`worstSeverity!(
                st;scopeVal;count scopeRows;breachCnt;warnCnt;maxUtil;worstSev)];
        rIdx+:1];
    resultRows
 };

.limitreport.summaryByMetric:{[limitCheckResult]
    metricNames:distinct limitCheckResult`metricName;
    resultRows:();
    mIdx:0;
    while[mIdx<count metricNames;
        mn:metricNames mIdx;
        metricMask:(limitCheckResult`metricName)=mn;
        metricRows:limitCheckResult where metricMask;
        sevCol:metricRows`severity;
        breachCnt:sum sevCol=`breach;
        warnCnt:sum sevCol=`warning;
        utilCol:metricRows`utilisation;
        validUtils:utilCol where not null utilCol;
        maxUtil:$[0<count validUtils;max validUtils;0Nf];
        worstSev:`OK;
        if[0<sum sevCol=`error; worstSev:`error];
        if[0<sum sevCol=`breach; worstSev:`breach];
        if[(worstSev=`OK) and 0<sum sevCol=`warning; worstSev:`warning];
        resultRows:resultRows,enlist `metricName`limitCount`breachCount`warningCount`maxUtilisation`worstSeverity!(
            mn;count metricRows;breachCnt;warnCnt;maxUtil;worstSev);
        mIdx+:1];
    resultRows
 };

.limitreport.limitDashboard:{[limitCheckResult]
    sevCol:limitCheckResult`severity;
    totalLimits:count limitCheckResult;
    okCnt:sum sevCol=`OK;
    warnCnt:sum sevCol=`warning;
    breachCnt:sum sevCol=`breach;
    errCnt:sum sevCol=`error;
    utilCol:limitCheckResult`utilisation;
    validUtils:utilCol where not null utilCol;
    maxUtil:$[0<count validUtils;max validUtils;0Nf];
    worstSev:`OK;
    if[errCnt>0; worstSev:`error];
    if[breachCnt>0; worstSev:`breach];
    if[(worstSev=`OK) and warnCnt>0; worstSev:`warning];
    `totalLimits`okCount`warningCount`breachCount`errorCount`worstSeverity`maxUtilisation`status`errorMessage!(
        totalLimits;okCnt;warnCnt;breachCnt;errCnt;worstSev;maxUtil;`OK;"")
 };

.limitreport.exportLimitReports:{[limitCheckResult;outputDirectory;reportLabel]
    / Export placeholder - returns report paths
    breaches:.limitreport.breachReport limitCheckResult;
    warnings:.limitreport.warningReport limitCheckResult;
    dashboard:.limitreport.limitDashboard limitCheckResult;
    `breachCount`warningCount`totalLimits`status!(count breaches;count warnings;count limitCheckResult;`OK)
 };
