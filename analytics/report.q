/ report.q - desk-style reporting: summaries, error filtering, CSV export

/ --- Public ---

.report.portfolioSummary:{[pricingResult]
    totalTrades:count pricingResult;
    okTrades:sum pricingResult[`status]=`OK;
    errorTrades:totalTrades-okTrades;
    okRows:pricingResult where pricingResult[`status]=`OK;
    totalNotional:0f;
    if[0<count okRows; totalNotional:sum okRows`notionalPrice];
    `totalTrades`okTrades`errorTrades`totalNotionalPrice!(
        totalTrades;okTrades;errorTrades;totalNotional)
 };

.report.riskSummary:{[greekResult]
    totalTrades:count greekResult;
    okRows:greekResult where greekResult[`status]=`OK;
    unsupportedRows:greekResult where greekResult[`status]=`UNSUPPORTED;
    okCount:count okRows;
    unsupportedCount:count unsupportedRows;
    portfolioDelta:0f; portfolioGamma:0f; portfolioVega:0f; portfolioRho:0f;
    if[okCount>0;
        portfolioDelta:sum okRows`delta;
        portfolioGamma:sum okRows`gamma;
        portfolioVega:sum okRows`vega;
        portfolioRho:sum okRows`rho];
    `totalTrades`okTrades`unsupportedTrades`portfolioDelta`portfolioGamma`portfolioVega`portfolioRho!(
        totalTrades;okCount;unsupportedCount;portfolioDelta;portfolioGamma;portfolioVega;portfolioRho)
 };

.report.scenarioSummary:{[scenarioResult]
    okRows:scenarioResult where scenarioResult[`status]=`OK;
    scenarioNames:distinct okRows`scenario;
    resultTable:();
    scenIdx:0;
    while[scenIdx<count scenarioNames;
        scenName:scenarioNames scenIdx;
        scRows:okRows where okRows[`scenario]=scenName;
        summaryRow:`scenario`tradeCount`totalNotionalPrice`totalNotionalPnL!(
            scenName;count scRows;sum scRows`notionalPrice;sum scRows`notionalPnL);
        resultTable:resultTable,enlist summaryRow;
        scenIdx+:1];
    resultTable
 };

.report.errorSummary:{[resultTable]
    resultTable where not resultTable[`status]=`OK
 };

.report.exportCsv:{[tbl;filePath]
    fileHandle:hsym `$filePath;
    fileHandle 0: csv 0: tbl;
    filePath
 };
