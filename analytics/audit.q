/ audit.q - audit records and error consolidation for batch runs

/ --- Public ---

.audit.createAuditRecord:{[runResult;valuationDate;runLabel]
    pricingResult:runResult`pricingResult;
    scenarioResult:runResult`scenarioResult;
    pnlExplainResult:runResult`pnlExplainResult;
    tradeCount:count pricingResult;
    pricedOkCount:sum pricingResult[`status]=`OK;
    errorCount:tradeCount-pricedOkCount;
    okPricingRows:pricingResult where pricingResult[`status]=`OK;
    totalPV:0f;
    if[0<count okPricingRows; totalPV:sum okPricingRows`notionalPrice];
    scenarioCount:count distinct scenarioResult`scenario;
    `valuationDate`runLabel`runTimestamp`tradeCount`pricedOkCount`errorCount`scenarioCount`totalPV`totalScenarioRows`totalPnlExplainRows!(
        valuationDate;
        runLabel;
        .z.p;
        tradeCount;
        pricedOkCount;
        errorCount;
        scenarioCount;
        totalPV;
        count scenarioResult;
        count pnlExplainResult)
 };

.audit.errorAudit:{[pricingResult;scenarioResult;pnlExplainResult]
    pricingErrors:pricingResult where not pricingResult[`status]=`OK;
    scenarioErrors:scenarioResult where not scenarioResult[`status]=`OK;
    pnlErrors:pnlExplainResult where not pnlExplainResult[`status]=`OK;
    `pricingErrorCount`scenarioErrorCount`pnlErrorCount`pricingErrors`scenarioErrors`pnlErrors!(
        count pricingErrors;count scenarioErrors;count pnlErrors;pricingErrors;scenarioErrors;pnlErrors)
 };
