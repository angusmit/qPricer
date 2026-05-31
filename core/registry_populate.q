/ core/registry_populate.q - register the EXISTING capabilities into the spine
/ (v0.68, Research OS R2). Loaded LAST (after every capability layer exists).
/ ----------------------------------------------------------------------------
/ Proves the spine is real + populated: each key existing capability is registered
/ through its per-kind registry with an accurate MANIFEST (declared I/O schema) and a
/ HANDLE to the UNCHANGED function. METADATA ONLY - no compute change, no canonical
/ number moves; the existing direct call paths (engine, backtest, exec) are untouched.
/ Every registration here must PASS .contracts.verify[] (if one does not, the contract
/ or the manifest is wrong - fix the manifest/contract, never the capability's compute).
/ Strategies remain in their own existing registry (.strategy.register) - unchanged.
/ ----------------------------------------------------------------------------

/ model / pricer: the engine's public FDM European-vanilla pricer.
.model.register[`blackScholesFdm; .engine.priceOption;
    `contractVersion`description`in`out!(
        1;
        "Crank-Nicolson / explicit FDM European-vanilla pricer (.engine.priceOption)";
        `spot`strike`expiry`rate`vol`notional!`float`float`float`float`float`float;
        `unitPrice`notionalPrice`method`modelName!`float`float`symbol`symbol)];

/ signal: the commodity signal path used by the backtests (carry / momentum / Kalman /
/ curve-residual, incl. the timeSeriesMomentum signal).
.signal.register[`commoditySignalPath; .strategy.path.commoditySignals;
    `contractVersion`description`in`out!(
        1;
        "commodity signal path from a forward-curve history (.strategy.path.commoditySignals)";
        `curveHistory`sigCfg!`table`dict;
        `path`kalmanParams`trainEndDate!`table`dict`date)];

/ calibrator: the two-factor Schwartz-Smith forward-curve calibration.
.calibrator.register[`schwartz2Curve; .commodity.calibrateCurve;
    `contractVersion`description`in`out!(
        1;
        "two-factor Schwartz-Smith forward-curve calibration (.commodity.calibrateCurve)";
        `marketCurve`modelType`calCfg!`table`symbol`dict;
        `calibratedParams`fitRmse`status!`dict`float`symbol)];

/ fillModel (execution): the generic daily fill-and-cost model.
.execution.fillModel.register[`dailyFillCost; .exec.fill;
    `contractVersion`description`in`out!(
        1;
        "generic daily fill-and-cost model: participation cap + slippage + costs (.exec.fill)";
        `order`ctx`cfg!`float`dict`dict;
        `filledQty`fillPrice`proportionalCost`slippageCost`fixedCost`totalCost!`float`float`float`float`float`float)];

-1 "registry_populate.q loaded - registered ",(string sum count each .registry.__store)," capabilities across ",(string count .registry.kinds[])," kinds";
