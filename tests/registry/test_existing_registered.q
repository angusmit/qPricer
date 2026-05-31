\l core/init.q
/ ============================================================================
/ test_existing_registered.q - the EXISTING capabilities registered at init load
/ (core/registry_populate.q) all conform to their contracts, and their fn handles
/ point at the real UNCHANGED functions. Research OS R2. No HDB read.
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};

/ --- each registered capability conforms to its kind's contract ---
chk[.model.conforms[`blackScholesFdm]; "the FDM pricer must conform to the model contract"];
chk[.signal.conforms[`commoditySignalPath]; "the commodity signal path must conform to the signal contract"];
chk[.calibrator.conforms[`schwartz2Curve]; "the schwartz2 curve calibrator must conform to the calibrator contract"];
chk[.execution.fillModel.conforms[`dailyFillCost]; "the daily fill model must conform to the fillModel contract"];

/ --- the spine is populated across all four kinds ---
chk[all `model`signal`calibrator`fillModel in .registry.kinds[]; "all four capability kinds must be registered"];

/ --- .contracts.verify[] passes for the real capabilities (no failures spine-wide) ---
v:.contracts.verify[];
chk[(count v)>=4; "verify must cover the >=4 registered capabilities"];
chk[all v`pass; "every registered capability must PASS .contracts.verify[]"];

/ --- the fn handles point at the real, unchanged functions (metadata only) ---
chk[(.model.get[`blackScholesFdm]`fn)~.engine.priceOption; "model handle must be .engine.priceOption"];
chk[(.signal.get[`commoditySignalPath]`fn)~.strategy.path.commoditySignals; "signal handle must be .strategy.path.commoditySignals"];
chk[(.calibrator.get[`schwartz2Curve]`fn)~.commodity.calibrateCurve; "calibrator handle must be .commodity.calibrateCurve"];
chk[(.execution.fillModel.get[`dailyFillCost]`fn)~.exec.fill; "fillModel handle must be .exec.fill"];

/ --- and the existing strategy registry is untouched (its own registry, not rebased) ---
chk[`gammaScalp in .strategy.registeredStrategies[]; "the strategy registry must still resolve gammaScalp"];

-1 "test_existing_registered: PASS";
