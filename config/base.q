/ config/base.q - baseline configuration values (v0.57, ARCHITECTURE.md section 2)
/ ----------------------------------------------------------------------------
/ The .cfg namespace of domain-keyed dicts. EVERY value here EXACTLY equals the
/ literal it replaced in the module it came from, so loading base alone (the
/ default / unset-QPRICER_ENV path) reproduces the prior hardcoded behaviour
/ byte-identically. config/{env}.q files (loaded after this by core/cfg.q when
/ QPRICER_ENV is set) may upsert a subset to override for dev / test / prod.
/ Types are preserved verbatim (a 200 stays a long, a 0f stays a float) so the
/ '~' dict comparisons in the suite are unchanged. Key ORDER matches the
/ original literals so the returned dicts compare equal.
/ ----------------------------------------------------------------------------

/ --- core layer: FDM pricing config (core/config.q) ---
/ defaultPricingConfig; withDefaults reads the trailing four keys from this dict.
.cfg.fdm:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

/ --- core layer: implied-volatility solver (core/config.q) ---
/ defaultImpliedVolConfig.
.cfg.iv:`lowerVolatilityBound`upperVolatilityBound`tolerance`maximumIterations!(
    0.0001;5.0;1e-8;100);

/ --- core layer: Monte Carlo defaults (core/montecarlo.q) ---
/ defaultMcConfig.
.cfg.mc:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(
    50000;50;42;1b;0b;0.95);

/ --- calibration layer: curve calibration (calibration/calibrateCurve.q) ---
/ defaultCalCfg: schwartz2 curve fit (vols/corr/rate fixed; kappa searched).
.cfg.calib.curve:`shortVolatility`longVolatility`correlation`riskFreeRate`meanReversionSpeedRange`gridSteps`refineRounds`refineShrink`kappaFloor!(
    0.30;0.15;0.30;0.02;0.1 3.0;25;5;0.3;0.05);

/ defaultSeriesCfg: per-date convenience-yield series (kappa fixed per date).
.cfg.calib.curveSeries:`kappa`shortVolatility`longVolatility`correlation`riskFreeRate!(
    1.0;0.30;0.15;0.30;0.02);

/ --- calibration layer: Kalman Schwartz-Smith (calibration/kalmanSchwartzSmith.q) ---
/ defaultParams: full param dict (lamChi=lamXi=0 fixed).
.cfg.calib.kalmanParams:`kappa`muXi`sigChi`sigXi`correlation`lamChi`lamXi`measSigma!(
    1.0;0.0;0.30;0.15;0.0;0.0;0.0;0.02);

/ defaultEstCfg: MLE bounds / init / refining-grid controls.
.cfg.calib.kalmanEst:`bounds`init`gridSteps`refineRounds`refineShrink`nSweeps!(
    `kappa`muXi`sigChi`sigXi`correlation`measSigma!((0.05 5.0);(-0.5 0.5);(0.01 1.5);(0.01 1.0);(-0.95 0.95);(0.0001 0.3));
    `kappa`muXi`sigChi`sigXi`correlation`measSigma!(1.0;0.0;0.30;0.15;0.0;0.02);
    11;4;0.4;5);

/ --- analytics layer: commodity model-report thresholds (analytics/commodityModelReport.q) ---
/ defaultGreekConfig: finite-difference bump sizes.
.cfg.analytics.greek:`spotBumpPct`volBump`meanReversionBump`longRunLogMeanBump`jumpIntensityBump`jumpMeanBump`correlationBump`useCentralDifference!(
    0.01;0.01;0.01;0.01;0.10;0.01;0.01;1b);

/ defaultDisagreementConfig: cross-model disagreement alert thresholds.
.cfg.analytics.disagreement:`priceRangeAbsThreshold`priceRangePctThreshold`primaryDeltaRangeThreshold`volatilityVegaRangeThreshold`scenarioPnlRangeThreshold`jumpSensitivityThreshold`minimumOkModels!(
    5f;0.10;10f;10f;1000f;5f;2);

/ defaultModelRiskLimitConfig: portfolio model-risk limit warning/breach levels.
.cfg.analytics.modelRiskLimit:`grossPriceRangeExposureWarning`grossPriceRangeExposureBreach`maxScenarioPnlRangeWarning`maxScenarioPnlRangeBreach`maxPrimarySensitivityRangeWarning`maxPrimarySensitivityRangeBreach`maxVolatilityVegaRangeWarning`maxVolatilityVegaRangeBreach`maxJumpSensitivityWarning`maxJumpSensitivityBreach`warningAlertCountWarning`warningAlertCountBreach`errorTradeCountWarning`errorTradeCountBreach!(
    5000f;10000f;1000f;2500f;25f;50f;20f;40f;5f;10f;5;10;1;3);

/ --- paths: default data + output directories ---
/ Default locations for real-data ingestion / report output. The library takes
/ directories as explicit arguments (apps supply them), so these are a documented
/ default registry for apps/pipeline; not read by tests, so they do not affect the
/ byte-identical guarantee. Repointing apps to read these is deferred.
.cfg.paths:`barchartRoot`crudeDir`gasDir`aaplOptionsDir`outputDir!(
    "data/barchart";
    "data/barchart/CRUDE";
    "data/barchart/GAS";
    "data/barchart/aapl/options_history";
    "output");
