/ var.q - VaR and Expected Shortfall risk engine (v0.27)

.var.validateConfidenceLevel:{[confidenceLevel]
    if[not confidenceLevel>0f; '"Confidence level must be > 0"];
    if[not confidenceLevel<1f; '"Confidence level must be < 1"];
 };

.var.valueAtRisk:{[pnlVector;confidenceLevel]
    .var.validateConfidenceLevel confidenceLevel;
    .riskdist.validatePnlVector pnlVector;
    sortedPnl:asc pnlVector;
    nObs:count sortedPnl;
    tailFraction:1f-confidenceLevel;
    tailIdx:floor tailFraction*nObs;
    tailIdx:0|tailIdx-1;
    neg sortedPnl tailIdx
 };

.var.expectedShortfall:{[pnlVector;confidenceLevel]
    .var.validateConfidenceLevel confidenceLevel;
    .riskdist.validatePnlVector pnlVector;
    sortedPnl:asc pnlVector;
    nObs:count sortedPnl;
    tailFraction:1f-confidenceLevel;
    tailCount:1|ceiling tailFraction*nObs;
    tailPnl:tailCount#sortedPnl;
    neg avg tailPnl
 };

.var.varExpectedShortfall:{[pnlVector;confidenceLevel]
    .var.validateConfidenceLevel confidenceLevel;
    .riskdist.validatePnlVector pnlVector;
    varVal:.var.valueAtRisk[pnlVector;confidenceLevel];
    esVal:.var.expectedShortfall[pnlVector;confidenceLevel];
    nObs:count pnlVector;
    tailFraction:1f-confidenceLevel;
    tailCount:1|ceiling tailFraction*nObs;
    `confidenceLevel`valueAtRisk`expectedShortfall`tailCount`observationCount`status`errorMessage!(
        confidenceLevel;varVal;esVal;tailCount;nObs;`OK;"")
 };

.var.varReport:{[pnlVector;confidenceLevels]
    resultRows:();
    clIdx:0;
    while[clIdx<count confidenceLevels;
          cl:confidenceLevels clIdx;
          rowResult:@[.var.varExpectedShortfall[pnlVector;];cl;{
              `confidenceLevel`valueAtRisk`expectedShortfall`tailCount`observationCount`status`errorMessage!(
                  0Nf;0Nf;0Nf;0N;0N;`ERROR;x)}];
          resultRows:resultRows,enlist rowResult;
          clIdx+:1];
    resultRows
 };

.var.varReportFromScenarioResult:{[scenarioResult;confidenceLevels]
    pnlDist:.riskdist.buildPnlDistribution scenarioResult;
    pnlVec:.riskdist.pnlVector pnlDist;
    .var.varReport[pnlVec;confidenceLevels]
 };

/ --- Risk contribution (standalone group VaR) ---

.var.riskContribution:{[pnlTable;groupColumn;confidenceLevel]
    .var.validateConfidenceLevel confidenceLevel;
    firstRow:pnlTable 0;
    tableCols:key firstRow;
    if[not groupColumn in tableCols; '"Missing group column: ",string groupColumn];
    if[not `scenarioName in tableCols; '"Missing scenarioName column"];
    if[not `pnl in tableCols; '"Missing pnl column"];
    groupNames:distinct pnlTable groupColumn;
    / Total VaR across all groups — rename pnl to scenarioPnl for compatibility
    totalScenario:();
    tIdx:0;
    while[tIdx<count pnlTable;
          rowData:pnlTable tIdx;
          totalScenario:totalScenario,enlist `scenarioName`scenarioPnl!(rowData`scenarioName;rowData`pnl);
          tIdx+:1];
    totalDist:.riskdist.buildPnlDistribution totalScenario;
    totalPnlVec:.riskdist.pnlVector totalDist;
    totalVarVal:.var.valueAtRisk[totalPnlVec;confidenceLevel];
    totalEsVal:.var.expectedShortfall[totalPnlVec;confidenceLevel];
    resultRows:();
    gIdx:0;
    while[gIdx<count groupNames;
          gn:groupNames gIdx;
          mask:(pnlTable groupColumn)=gn;
          groupRows:pnlTable where mask;
          / Rename pnl to scenarioPnl for buildPnlDistribution compatibility
          groupScenario:();
          rIdx:0;
          while[rIdx<count groupRows;
                rowData:groupRows rIdx;
                groupScenario:groupScenario,enlist `scenarioName`scenarioPnl!(rowData`scenarioName;rowData`pnl);
                rIdx+:1];
          groupDist:.riskdist.buildPnlDistribution groupScenario;
          groupPnlVec:.riskdist.pnlVector groupDist;
          groupVarVal:.var.valueAtRisk[groupPnlVec;confidenceLevel];
          groupEsVal:.var.expectedShortfall[groupPnlVec;confidenceLevel];
          varContribPct:$[totalVarVal>0f;groupVarVal%totalVarVal;0Nf];
          esContribPct:$[totalEsVal>0f;groupEsVal%totalEsVal;0Nf];
          resultRows:resultRows,enlist `groupName`confidenceLevel`groupVaR`groupExpectedShortfall`totalVaR`totalExpectedShortfall`varContributionPct`esContributionPct`scenarioCount`status`errorMessage!(
              gn;confidenceLevel;groupVarVal;groupEsVal;totalVarVal;totalEsVal;varContribPct;esContribPct;
              count groupDist;`OK;"");
          gIdx+:1];
    resultRows
 };

.var.componentVarReport:{[pnlTable;groupColumn;confidenceLevel]
    .var.riskContribution[pnlTable;groupColumn;confidenceLevel]
 };
