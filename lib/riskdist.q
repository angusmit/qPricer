/ riskdist.q - PnL distribution utilities (v0.27)

.riskdist.validatePnlVector:{[pnlVector]
    if[0=count pnlVector; '"PnL vector is empty"];
    if[any null pnlVector; '"PnL vector contains nulls"];
 };

.riskdist.validateScenarioResult:{[scenarioResult]
    if[0=count scenarioResult; '"Scenario result is empty"];
    firstRow:scenarioResult 0;
    tableCols:key firstRow;
    if[not `scenarioName in tableCols; '"Missing scenarioName column"];
    if[not `scenarioPnl in tableCols; '"Missing scenarioPnl column"];
 };

.riskdist.buildPnlDistribution:{[scenarioResult]
    .riskdist.validateScenarioResult scenarioResult;
    scenarioNames:distinct scenarioResult`scenarioName;
    resultRows:();
    sIdx:0;
    while[sIdx<count scenarioNames;
        sn:scenarioNames sIdx;
        mask:(scenarioResult`scenarioName)=sn;
        matchRows:scenarioResult where mask;
        totalPnl:sum matchRows`scenarioPnl;
        resultRows:resultRows,enlist `scenarioName`pnl`status`errorMessage!(sn;totalPnl;`OK;"");
        sIdx+:1];
    resultRows
 };

.riskdist.aggregatePnlByScenario:{[scenarioResult]
    .riskdist.buildPnlDistribution scenarioResult
 };

.riskdist.pnlVector:{[pnlDistribution]
    pnlDistribution`pnl
 };

.riskdist.sortPnl:{[pnlVector] asc pnlVector};

.riskdist.lossVector:{[pnlVector] neg pnlVector};

.riskdist.distributionSummary:{[pnlVector]
    .riskdist.validatePnlVector pnlVector;
    sortedPnl:asc pnlVector;
    nObs:count pnlVector;
    medianIdx:floor nObs%2;
    losses:pnlVector where pnlVector<0f;
    lossMeanVal:$[0<count losses;neg avg losses;0f];
    lossMaxVal:$[0<count losses;neg min losses;0f];
    `observationCount`meanPnl`medianPnl`minPnl`maxPnl`pnlStd`lossMean`lossMax`status`errorMessage!(
        nObs;avg pnlVector;sortedPnl medianIdx;min pnlVector;max pnlVector;dev pnlVector;
        lossMeanVal;lossMaxVal;`OK;"")
 };

.riskdist.worstLossScenarios:{[pnlDistribution;topCountVal]
    pnlCol:pnlDistribution`pnl;
    sortIdx:iasc pnlCol;
    takeCount:topCountVal&count pnlDistribution;
    resultRows:();
    rankIdx:0;
    while[rankIdx<takeCount;
        origIdx:sortIdx rankIdx;
        rowData:pnlDistribution origIdx;
        resultRows:resultRows,enlist `scenarioName`pnl`loss`rank`status`errorMessage!(
            rowData`scenarioName;rowData`pnl;neg rowData`pnl;rankIdx+1;`OK;"");
        rankIdx+:1];
    resultRows
 };
