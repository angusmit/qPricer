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
.cfg.paths:`barchartRoot`crudeDir`gasDir`aaplOptionsDir`outputDir`hdb!(
    "data/barchart";
    "data/barchart/CRUDE";
    "data/barchart/GAS";
    "data/barchart/aapl/options_history";
    "output";
    "data/hdb");

/ --- execution layer: daily fill-and-cost simulation (execution/execution.q) ---
/ Realism is OFF by default so the wired backtest is byte-identical to the legacy
/ flat-cost path: proportional cost only (at the strategy's own txnCostRate, which
/ .exec.__resolve injects over proportionalRate), zero slippage, zero fixed cost,
/ no participation cap (0n -> full fill). slippageBps / fixedPerTrade /
/ participationCap / impactCoef / volScaledSlippage are the opt-in realism levers.
/ Types matter for byte-identity: proportionalRate/slippageBps/impactCoef/
/ fixedPerTrade/participationCap are floats, volScaledSlippage a bool, costMode a sym.
.cfg.exec:`costMode`proportionalRate`slippageBps`impactCoef`volScaledSlippage`fixedPerTrade`participationCap!(
    `proportional;0.0005;0f;0f;0b;0f;0n);

/ --- portfolio layer: strategy allocator (portfolio/portfolio.q) ---
/ Defaults for .alloc.* : long-only + fully-invested (a strategy-sleeve allocator),
/ weightCap off (1f), no turnover penalty, riskAversion 1 (meanVariance), no
/ covariance shrinkage, 300 ERC iterations, 252-day annualisation. riskParity is
/ the recommended method (covariance-only, no return forecast); equalWeight is the baseline.
.cfg.alloc:`longOnly`fullyInvested`weightCap`turnoverPenalty`riskAversion`shrinkage`rpMaxIter`annualizationDays!(
    1b;1b;1f;0f;1f;0f;300;252f);

/ --- regime layer: Market State Engine thresholds (regime/regime.q) ---
/ deferredIdx = which tenor is the 'deferred' leg for the curve slope (5 -> the 6th
/ contract, ~CL6). flatSlopeThreshold: |relative slope| below this is `flat. volLookback:
/ realized-vol window (days). pctLookback: rolling-percentile window. rollNearDays: front
/ days-to-expiry <= this is `near. low/high + thin/deep percentile cutoffs for vol/liquidity.
/ seasonPos/NegThreshold + seasonCfg drive the seasonal-phase bucket via signals/seasonality.
.cfg.regime:`deferredIdx`flatSlopeThreshold`volLookback`pctLookback`rollNearDays`lowPctThreshold`highPctThreshold`thinPctThreshold`deepPctThreshold`seasonPosThreshold`seasonNegThreshold`seasonCfg`annualizationDays!(
    5;0.005;21;252;10;0.33;0.67;0.33;0.67;0.01;-0.01;`amplitude`phaseYears!(0.03;0f);252f);

/ --- governance layer: research-governance gates (gov/gov.q) ---
/ Research OS R3. costHurdleSharpe: Gate 1 net-of-execution annualised-Sharpe floor.
/ dsrThreshold: Gate 2 deflated-Sharpe pass bar (Bailey & Lopez de Prado standard 0.95).
/ validEdgeSources: the edgeSource enum Gate 0 accepts. exec: a REALISTIC .exec cost
/ sub-config (slippageBps etc.) the wrapper/demo threads into the backtest so the
/ cascade judges NET returns (NOT the frictionless default). walkForward: Gate 3 causal
/ split scheme (reuses .strategy.commodityBT.__splits) + the OOS sign-stability floor.
.cfg.gov:`costHurdleSharpe`dsrThreshold`validEdgeSources`annualizationDays`exec`walkForward!(
    0.25;
    0.95;
    `riskPremium`mispricing`info;
    252f;
    `slippageBps`fixedPerTrade!(2f;0f);
    `scheme`trainSpan`testSpan`maxSplits`minFolds`oosSharpeFloor`minPassFraction!(`rolling;126;63;6;3;0f;0.6));

/ --- backtest layer: per-strategy default configs (backtest/strategy.q) ---
/ Each .strategy.<name>.defaultConfig returns its dict below. Values, TYPES and
/ key order are verbatim from the prior function literals (1f%252f kept as the
/ exact expression; longs like rebalanceInterval=1 stay long; 0n / 0Nd / ()!()
/ preserved), so behaviour is byte-identical.

.cfg.strategy.gammaScalp:`optionSide`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears`useModelDelta!(
    `long;`interval;1;0.05;0f;0f;1f%252f;1b);

.cfg.strategy.shortVariance:`forecastVol`entryMargin`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears!(
    0.20;0.02;`interval;1;0.05;0f;0f;1f%252f);

.cfg.strategy.calendarRoll:`spreadType`optionType`frontTenorYears`backTenorYears`rollThresholdYears`rollBackLeg`hedgeDelta`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears!(
    `longCalendar;`call;0.05;0.20;0.005;1b;1b;`interval;1;0.05;0f;0f;1f%252f);

.cfg.strategy.riskReversal:`riskReversalDirection`wingOffsetPct`skewSlope`fairSkew`skewMargin`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears`hedgeDelta!(
    `auto;0.05;-0.5;-0.2;0.05;`interval;1;0.05;0f;0f;1f%252f;1b);

.cfg.strategy.modelDisagreement:`modelAVolBump`modelBVolBump`disagreementThreshold`direction`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears`hedgeDelta!(
    0f;0.05;0.30;`auto;`interval;1;0.05;0f;0f;1f%252f;1b);

.cfg.strategy.deltaVegaHedge:`hedgeOptionType`hedgeStrikeRatio`hedgeOptionExpiryYears`vegaRebalanceMode`vegaRebalanceInterval`vegaBand`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears!(
    `call;1.05;0n;`interval;1;0.05;`interval;1;0.05;0f;0f;1f%252f);

.cfg.strategy.longVol:`forecastVol`entryMargin`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears!(
    0.30;0.02;`interval;1;0.05;0f;0f;1f%252f);

.cfg.strategy.collarTailHedge:`mode`putStrikePct`callStrikePct`premiumBudgetPct`txnCostRate`financingRate`stepYears!(
    `collar;0.05;0.05;0.005;0f;0f;1f%252f);

.cfg.strategy.putRatioBackspread:`shortStrikePct`longStrikePct`ratioN`hedgeDelta`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears!(
    0f;0.05;2f;1b;`interval;1;0.05;0f;0f;1f%252f);

.cfg.strategy.ironCondor:`shortPutPct`longPutPct`shortCallPct`longCallPct`hedgeDelta`entryVolThreshold`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears!(
    0.05;0.10;0.05;0.10;0b;0f;`interval;1;0.05;0f;0f;1f%252f);

.cfg.strategy.barrierHedge:`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears!(
    `interval;1;0.05;0f;0f;1f%252f);

.cfg.strategy.jumpPremium:`jumpModelName`jumpModelParams`thresholdType`premiumThreshold`direction`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears!(
    `merton;
    `volatility`jumpIntensity`jumpMean`jumpVolatility`riskFreeRate`dividendYield!(0.20;0.5;-0.05;0.20;0.02;0f);
    `absolute;0.5;`auto;`interval;1;0.05;0f;0f;1f%252f);

.cfg.strategy.dispersion:`pathBundle`indexVol`direction`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears!(
    ()!();0.20;`shortIndex;`interval;1;0.05;0f;0f;1f%252f);

.cfg.strategy.powerSpikeCapture:`curveBundle`callStrikePct`callExpiry`callVol`deviationThreshold`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears!(
    ()!();0.05;0.10;0.40;2f;`interval;1;0.05;0f;0f;1f%252f);

.cfg.strategy.commodityCalendar:`curveBundle`nearTenor`farTenor`rollTriggerTenor`txnCostRate`financingRate`stepYears!(
    ()!();0.10;0.30;0.05;0f;0f;1f%252f);

.cfg.strategy.sparkSpread:`curveBundle`leg1Name`leg2Name`heatRate`strike`optType`vol1`vol2`correlation`expiry`hedgeEnabled`bumpSize`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears!(
    ()!();`power;`gas;8f;0f;`call;0.45;0.35;0.4;0.25;1b;0.01;`interval;1;0.05;0f;0f;1f%252f);

.cfg.strategy.crackSpread:`curveBundle`leg1Name`leg2Name`crackRatio`strike`optType`vol1`vol2`correlation`expiry`hedgeEnabled`bumpSize`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears!(
    ()!();`product;`crude;1f;0f;`call;0.35;0.30;0.6;0.25;1b;0.01;`interval;1;0.05;0f;0f;1f%252f);

/ --- backtest layer: commodity strategy default configs (backtest/commodityStrategies.q) ---

.cfg.strategy.convenienceYieldCarry:`signalSource`carryMargin`riskFreeRate`allowShort`txnCostRate`annualizationDays`notional!(
    `convenienceYield;0.0;0.02;1b;0.0005;252f;1f);

.cfg.strategy.chiReversion:`entryZ`exitZ`txnCostRate`annualizationDays`notional!(
    1.0;0.3;0.0005;252f;1f);

.cfg.strategy.timeSeriesMomentum:`momentumMargin`txnCostRate`annualizationDays`notional!(
    0f;0.0005;252f;1f);

.cfg.strategy.twoTimescale:`revertWeight`trendWeight`entryZ`txnCostRate`annualizationDays`notional!(
    0.5;0.5;1.0;0.0005;252f;1f);

.cfg.strategy.storageCashCarry:`storageCostRate`carryMargin`txnCostRate`annualizationDays`notional`returnColumn`volScaleColumn!(
    0.01;0.0f;0.0005;252f;1f;`nearFarSpreadReturn;`volTargetScaleNearFar);

.cfg.strategy.carryMomentumCombo:`carryWeight`momentumWeight`riskFreeRate`carryMargin`momentumMargin`txnCostRate`annualizationDays`notional!(
    0.5;0.5;0.02;0.0f;0.0f;0.0005;252f;1f);

.cfg.strategy.curveRelativeValue:`minGap`txnCostRate`annualizationDays`notional`returnColumn`volScaleColumn!(
    0.0f;0.0005;252f;1f;`rvSpreadReturn;`volTargetScaleRV);

/ --- backtest layer: commodity signal-path config (backtest/commodityStrategies.q) ---
/ .strategy.path.__commodityDefaultSigCfg returns this dict.
.cfg.strategy.commoditySignals:`rollDaysBeforeExpiry`trainFraction`trainEndDate`momentumLookback`txnCostRate`targetVol`riskFreeRate`storageCostRate`annualizationDays`kalmanEstCfg`carryMargin`productTag`deseasonalize!(
    5;0.6;0Nd;20;0.0005;0.15;0.02;0.01;252f;`gridSteps`refineRounds`nSweeps!(7;3;3);0.0;`CRUDE;0b);
