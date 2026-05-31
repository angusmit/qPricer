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
/ R4 adds `analogue (the analogue-distance weights: weighted Euclidean over the three
/ percentile axes + a per-mismatch categorical penalty over the discrete states) and
/ `episodes (below). The analogue distance/ranking are pure; the episode spec is curated.
.cfg.regime:`deferredIdx`flatSlopeThreshold`volLookback`pctLookback`rollNearDays`lowPctThreshold`highPctThreshold`thinPctThreshold`deepPctThreshold`seasonPosThreshold`seasonNegThreshold`seasonCfg`annualizationDays`analogue!(
    5;0.005;21;252;10;0.33;0.67;0.33;0.67;0.01;-0.01;`amplitude`phaseYears!(0.03;0f);252f;
    `wSlope`wVol`wVolume`catPenalty!(1f;1f;1f;0.5));

/ Curated regime-episode spec (Research OS R4). ONLY in-HDB-coverage crude windows
/ (~2018.12 -> 2026) - out-of-data lessons (2008, 2014-16) live in docs/REGIME_LIBRARY.md
/ as narrative, never as fabricated fingerprints. driversKey -> the doc section; riskMemory
/ is a one-line summary (the full failure-mode writeup is in REGIME_LIBRARY.md).
/ scripts/build_regime_library.q computes each episode's dominant fingerprint from `regimes`.
.cfg.regime[`episodes]:([]
    episodeId:`crude_2020_covid`crude_2022_energyShock`crude_2023_rangebound;
    commodity:`CRUDE`CRUDE`CRUDE;
    dateFrom:2020.02.20 2022.02.24 2023.05.01;
    dateTo:2020.06.30 2022.08.31 2023.12.31;
    label:`crude_2020_covid`crude_2022_energyShock`crude_2023_rangebound;
    driversKey:`covid2020`energyShock2022`opecCuts2023;
    riskMemory:(
        "demand collapse + storage saturation: long flat-price and short-vol blew up, WTI printed negative 2020-04-20, calendar spreads dislocated to super-contango, liquidity vanished";
        "Russia/Ukraine supply shock: sharp backwardation + a vol spike, short-vol and mean-reversion shorts ran over, late momentum chasers round-tripped";
        "range-bound chop under OPEC+ cuts: trend/momentum whipsawed, breakouts faded, only modest carry paid"));

/ --- governance layer: research-governance gates (gov/gov.q) ---
/ Research OS R3. costHurdleSharpe: Gate 1 net-of-execution annualised-Sharpe floor.
/ dsrThreshold: Gate 2 deflated-Sharpe pass bar (Bailey & Lopez de Prado standard 0.95).
/ validEdgeSources: the edgeSource enum Gate 0 accepts. exec: a REALISTIC .exec cost
/ sub-config (slippageBps etc.) the wrapper/demo threads into the backtest so the
/ cascade judges NET returns (NOT the frictionless default). walkForward: Gate 3 causal
/ split scheme (reuses .strategy.commodityBT.__splits) + the OOS sign-stability floor.
/ R3b adds zones + holdoutSharpe: zones split each commodity's sorted dates into
/ train/validate/holdout (holdout = the most-recent 1-trainFrac-validateFrac slice,
/ out-of-sample in TIME); gates 0-3 see train+validate only. holdoutSharpe is the Gate 4
/ net-annualised-Sharpe bar on the sealed holdout (one-shot).
.cfg.gov:`costHurdleSharpe`dsrThreshold`validEdgeSources`annualizationDays`exec`walkForward`zones`holdoutSharpe!(
    0.25;
    0.95;
    `riskPremium`mispricing`info;
    252f;
    `slippageBps`fixedPerTrade!(2f;0f);
    `scheme`trainSpan`testSpan`maxSplits`minFolds`oosSharpeFloor`minPassFraction!(`rolling;126;63;6;3;0f;0.6);
    `trainFrac`validateFrac!(0.6;0.2);
    0.25);
/ v0.71 connective wiring: how many nearest historical episodes the regime risk-memory
/ skeptic surfaces on each verdict (informational annotation; changes no gate pass/fail).
.cfg.gov[`skepticN]:2;

/ --- templates layer: relativeValue problem template (templates/relative_value.q), R6 ---
/ lookback/entryZ drive the causal mean-reversion z-score (fade extremes); notional scales
/ the position; txnCostRate feeds the fill model; deferredIdx picks the deferred leg for the
/ calendar spread; maxAr1 is the spread-stationarity bar (AR(1) coef below this = mean-reverting).
.cfg.templates.rv:`lookback`entryZ`notional`txnCostRate`deferredIdx`maxAr1!(60;1.5;1f;0.0005;1;0.98);

/ --- factor layer: curve PCA (factor/factor.q), Research OS R8 ---
/ k = number of principal components (level/slope/curvature); tol/maxIter pin the deterministic
/ power iteration; nMaturities = the first-M tenor-ranked contracts forming the curve panel;
/ minLoadingStability = the factor-stability gate bar (cosine of first- vs second-half top loading);
/ residualMaturity = the 0-based maturity whose cumulative residual the factorRelativeValue template trades.
/ Do NOT tune k or the lookback to pass the gates - that is exactly where factor PCA overfits.
.cfg.factor:`k`tol`maxIter`nMaturities`minLoadingStability`residualMaturity!(3;1e-10;1000;5;0.9;2);

/ --- state layer: as-of market-state universe/liquidity filters (state/state.q), R9 ---
/ expiryStrictlyAfter: a contract is tradable as-of asOf iff its expiry is strictly after asOf
/ (1b) or on/after (0b). minVolume: a declared liquidity-filter slot (0 = no filter; the curve
/ engine R10 + roll R11 will use it). The as-of accessor itself takes only date<=asOf rows.
.cfg.state:`expiryStrictlyAfter`minVolume!(1b;0f);

/ --- curve layer: the curve engine (curve/curve.q), Research OS R10 ---
/ deferredIdx + flatThreshold are SOURCED FROM .cfg.regime so the curve engine's slope +
/ backwardation/contango/flat classification cannot drift from regime/'s (the agreement
/ precondition for a later rebase). minVolume = the liquidity floor for the tenor filter
/ (0 = keep all positive-price tenors).
.cfg.curve:`deferredIdx`flatThreshold`minVolume!(.cfg.regime`deferredIdx; .cfg.regime`flatSlopeThreshold; 0f);

/ --- season layer: curve/spread seasonality (season/season.q), Research OS R14 ---
/ deferredIdx = the deferred leg of the seasonal front-deferred spread (2 -> M1-M3, the calendar
/ spread R16 fades). minObs = the minimum SAME-CALENDAR-MONTH historical observations required for a
/ z-score; below it the z is null + lowConfidence (thin same-month history early in the data is noisy).
/ The seasonal stats are CAUSAL (same-calendar-month up to asOf, through R9's door) - never full-sample.
.cfg.season:`deferredIdx`minObs!(2;3);

/ --- carry layer: carry / storage economics (carry/carry.q), Research OS R14 ---
/ riskFreeRate (r) + storageCost are ASSUMED (config, NOT market-observed) - the convenience-yield +
/ cash-and-carry outputs are assumption-dependent and flagged as such. The implied carry comes from the
/ curve (we have it); the inventory-tightness is a PROXY from the degree of backwardation (no real
/ inventory data). carry reuses .cfg.curve.deferredIdx (the CL1-CLk leg) so impliedCarry == -rollYield.
.cfg.carry:`riskFreeRate`storageCost!(0.04;0.02);

/ --- roll layer: per-commodity roll rules (roll/roll.q), Research OS R11 ---
/ One rule per commodity (+ a `default fallback): type (days_before_expiry [the robust default,
/ needs expiry only] / volume_switch [as-of volume] / fixed_calendar [rollDays-th day of the
/ expiry month] / oi_switch [recognised but UNAVAILABLE - no OI column]), rollDays (days before
/ expiry to roll, or the calendar day), W (the roll-window length, in days, over which old->new
/ blends + the trailing-volume window), priceSource, and the continuous back-adjust method.
.cfg.rolls:`CRUDE`GAS`default!(
    `type`rollDays`W`priceSource`method!(`days_before_expiry;10;5;`settle;`difference);
    `type`rollDays`W`priceSource`method!(`days_before_expiry;14;5;`settle;`difference);
    `type`rollDays`W`priceSource`method!(`days_before_expiry;10;5;`settle;`difference));

/ --- replay layer: the event-driven replay engine (backtest/replay.q), Research OS R12 ---
/ The realism levers for replay mode, ALL OFF / ZERO by default so a frictionless replay uses the
/ UNCHANGED .exec.fill and reproduces research-mode PnL to tolerance (byte-identical canonicals).
/ participationRate = cap a daily fill at this fraction of the as-of bar volume (<=0 -> no cap);
/ remainderPolicy = what to do with the unfilled remainder under a cap (`carry to next day / `drop);
/ rollPenaltyBps = extra adverse slippage (bps of traded notional) on a roll date (R11 flags it);
/ financingRate = annualised carry on the held position value (0 -> none). rollW = the R11 roll
/ window for the replay's active-contract rule (1 -> switch exactly at rollDays, tracking the
/ research-mode heldSeq so the frictionless replay is faithful).
.cfg.replay:`participationRate`remainderPolicy`rollPenaltyBps`financingRate`annualizationDays`rollW!(
    0f;`drop;0f;0f;252f;1);

/ --- evidence layer: the deterministic evidence audit (evidence/evidence.q), Research OS R13 ---
/ reconcileTol = the absolute tolerance for the PnL-ties + book-ties reconciliation checks (the
/ replay folds floating-point sequentially, so allow a few ULPs). requireCosts = require a positive
/ transaction cost on a run that has fills (set 0b only for a genuinely zero-cost / frictionless run,
/ so the costs-applied check does not false-fail it). The audit BITES: any check fail -> overall fail.
.cfg.evidence:`reconcileTol`requireCosts!(1e-8; 1b);

/ --- attribution layer: PnL explain + bucketed curve risk (attribution/attribution.q), R15 ---
/ reconcileTol = the tolerance for the decomposition reconciliation (level+slope+curvature+carry+
/ residual == the run's realized total; the residual is the PLUG that absorbs the unexplained curve
/ move + costs). buckets = the tenor bucket edges (years) for the bucketed curve delta.
.cfg.attribution:`reconcileTol`buckets!(1e-8; 0 0.25 0.5 1 2 5f);

/ --- cards layer: curated model cards (cards/cards.q), Research OS R5 ---
/ One structured card per KEY capability (the 4 R2-registered capabilities + the key
/ strategies). capabilityName MATCHES the R2/strategy registry name. edgeSource is one of
/ riskPremium/structural/informational for signals/strategies, `na for pricers/fill models.
/ govHypoId links the card to its gov-ledger record (` if not gated). riskMemoryKey links
/ the regime library (`na if none). The full narrative is in docs/MODEL_CARDS.md.
.cfg.cards:([]
    cardId:`card_blackScholesFdm`card_commoditySignalPath`card_schwartz2Curve`card_dailyFillCost`card_gammaScalp`card_shortVariance`card_timeSeriesMomentum;
    capabilityKind:`model`signal`calibrator`fillModel`strategy`strategy`strategy;
    capabilityName:`blackScholesFdm`commoditySignalPath`schwartz2Curve`dailyFillCost`gammaScalp`shortVariance`timeSeriesMomentum;
    version:`v1`v1`v1`v1`v1`v1`v1;
    intendedUse:(
        "Price European-vanilla options via Crank-Nicolson / explicit FDM (the engine's reference pricer)";
        "Generate carry / momentum / Kalman / curve-residual signals from a forward-curve history";
        "Calibrate a two-factor Schwartz-Smith forward curve to a market curve";
        "Turn a target order into a simulated fill + explicit costs (participation cap + slippage)";
        "Long-gamma delta-hedged option position harvesting realized-over-implied vol";
        "Sell a delta-hedged ATM straddle to harvest the variance risk premium when IV is rich";
        "Long/short the front by the sign of trailing curve momentum (trend following)");
    assumptions:(
        "lognormal diffusion; constant vol/rate over the grid; European exercise; no early-exercise/barrier+localvol combos";
        ">=3 curve dates; curve history columns asofDate/tenor/price/contractYM/expiry; Kalman estimated train-only";
        ">=3 positive-price curve points; lognormal domain; schwartz2 model family";
        "daily settle prices only (no intraday/LOB); cost = proportional + slippage(bps of notional) + fixed";
        "continuous delta hedging at the configured band; realized vol exceeds implied; liquid underlying";
        "implied richer than the realized-vol forecast; straddle held to expiry; hedge slippage modest";
        "trend persistence in the front; vol-targeted sizing; turnover survives realistic execution cost");
    edgeSource:`na`riskPremium`na`na`riskPremium`riskPremium`riskPremium;
    regimeApplicability:(
        "any (a pricer, regime-agnostic)";
        "trending regimes (backwardation / strong directional); weak in range-bound chop";
        "any liquid forward curve (calibration utility, regime-agnostic)";
        "any (an execution model, regime-agnostic)";
        "high-realized-vol regimes; loses theta in calm/range-bound markets";
        "calm / mean-reverting vol regimes; blows up in vol spikes (2020/2022 type)";
        "supply-driven backwardation trends; whipsaws in OPEC-cut-defended ranges");
    riskMemoryKey:`na`energyShock2022`na`na`na`na`energyShock2022;
    govHypoId:(`;`mom_crude_ts;`;`;`;`;`mom_crude_ts);
    owner:(
        "qFDM core";"qFDM commodity desk";"qFDM commodity desk";"qFDM execution";
        "qFDM vol desk";"qFDM vol desk";"qFDM commodity desk");
    asOf:7#2026.05.31);

/ Research OS R8: cards for the first capability/template added on the spine. curvePCA is a
/ `factor capability (edge `na); factorRelativeValue is the template (edge `structural), carded
/ WITH a populated failure-mode link (riskMemoryKey covid2020 - the regime where the curve's
/ factor structure most violently broke) so it is gateable via .cards.gatedRun.
.cfg.cards:.cfg.cards upsert `cardId`capabilityKind`capabilityName`version`intendedUse`assumptions`edgeSource`regimeApplicability`riskMemoryKey`govHypoId`owner`asOf!(
    `card_curvePCA;`factor;`curvePCA;`v1;
    "PCA of the curve-change panel into level/slope/curvature + residuals (deterministic power iteration)";
    "k factors capture the structure; loadings stable over the window; deterministic (fixed init/tol/sign)";
    `na;
    "any liquid forward curve (an analytical capability, regime-agnostic)";
    `na;`;"qFDM quant research";2026.05.31);
.cfg.cards:.cfg.cards upsert `cardId`capabilityKind`capabilityName`version`intendedUse`assumptions`edgeSource`regimeApplicability`riskMemoryKey`govHypoId`owner`asOf!(
    `card_factorRelativeValue;`template;`factorRelativeValue;`v1;
    "Fade the curve's cumulative residual from its k-factor PCA shape (factor-structure reversion)";
    "factor structure stable across the window; the residual mean-reverts (not microstructure); k + lookback NOT tuned to the gates";
    `structural;
    "stable-factor-structure regimes; breaks down when the curve dislocates (storage saturation / supply shock)";
    `covid2020;`rv_factor_crude;"qFDM quant research";2026.05.31);
/ Research OS R9: the as-of accessor + Market State (a data-access capability, edge `na).
.cfg.cards:.cfg.cards upsert `cardId`capabilityKind`capabilityName`version`intendedUse`assumptions`edgeSource`regimeApplicability`riskMemoryKey`govHypoId`owner`asOf!(
    `card_marketState;`state;`marketState;`v1;
    "Point-in-time Market State for (asOf, commodity): the single as-of door to history + the state object";
    "the HDB expiry field is correct; futures' simple date<=asOf path (no revisable data); the accessor is the only door (existing direct readers not yet rebased)";
    `na;
    "any (a point-in-time data accessor, regime-agnostic)";
    `na;`;"qFDM evidence layer";2026.05.31);
/ Research OS R10: the curve engine (a derived-curve capability, edge `na).
.cfg.cards:.cfg.cards upsert `cardId`capabilityKind`capabilityName`version`intendedUse`assumptions`edgeSource`regimeApplicability`riskMemoryKey`govHypoId`owner`asOf!(
    `card_curveEngine;`curve;`curveEngine;`v1;
    "The curve engine: clean as-of curve + spreads + derived features (roll yield / slope / curvature / classification) + shocks + immutable snapshots";
    "reads through R9's as-of door (point-in-time by construction); slope+classification match regime/'s convention; snapshots assume a deterministic curve";
    `na;
    "any (a derived-curve capability, regime-agnostic)";
    `na;`;"qFDM evidence layer";2026.05.31);
/ Research OS R11: roll discipline (an as-of roll-mapping capability, edge `na).
.cfg.cards:.cfg.cards upsert `cardId`capabilityKind`capabilityName`version`intendedUse`assumptions`edgeSource`regimeApplicability`riskMemoryKey`govHypoId`owner`asOf!(
    `card_rollEngine;`roll;`rollEngine;`v1;
    "As-of-only active-contract roll rules + roll events + a (clearly-labeled, analytics-only) continuous back-adjusted view";
    "roll decisions use ONLY as-of data (through R9's door); trade the active contracts, NEVER the continuous series (non-point-in-time); no open interest in the HDB";
    `na;
    "any (an as-of roll-mapping capability, regime-agnostic)";
    `na;`;"qFDM evidence layer";2026.05.31);
/ Research OS R14: seasonality + carry FEATURE capabilities (the signals R16 is built on).
.cfg.cards:.cfg.cards upsert `cardId`capabilityKind`capabilityName`version`intendedUse`assumptions`edgeSource`regimeApplicability`riskMemoryKey`govHypoId`owner`asOf!(
    `card_curveSeasonality;`season;`curveSeasonality;`v1;
    "Causal curve/spread seasonality: same-calendar-month + same-contract-month z, seasonal factor, deseasonalised level, seasonal slope - the signal R16's calendar-spread mean-reversion fades";
    "the same-month stat is CAUSAL (same-calendar-month up to asOf, through R9's door) - never full-sample; thin same-month history early in the data is noisy (mitigated by min-N -> null/low-confidence); seasonal patterns break in structural shifts; the causal trailing stats use less data than full-sample";
    `structural;
    "commodities with a genuine seasonal curve (gas/power strongly, crude weakly); breaks down in structural shifts";
    `na;`;"qFDM evidence layer";2026.05.31);
.cfg.cards:.cfg.cards upsert `cardId`capabilityKind`capabilityName`version`intendedUse`assumptions`edgeSource`regimeApplicability`riskMemoryKey`govHypoId`owner`asOf!(
    `card_carryEconomics;`carry;`carryEconomics;`v1;
    "Carry/storage economics from the curve: implied carry, convenience yield + cash-and-carry fair value (assumption-dependent), carry signal, inventory-tightness proxy";
    "implied carry comes from the curve; convenience yield + cash-and-carry fair value depend on the ASSUMED r+storage (.cfg.carry, NOT market-observed) and are flagged; the inventory-tightness is a PROXY from the degree of backwardation (no real inventory data); cash-and-carry assumes storability";
    `riskPremium;
    "storable commodities (the convenience-yield / cash-and-carry interpretation assumes storability)";
    `na;`;"qFDM evidence layer";2026.05.31);
/ Research OS R15: PnL attribution (a decomposition capability, edge `na - it explains a run, not trades).
.cfg.cards:.cfg.cards upsert `cardId`capabilityKind`capabilityName`version`intendedUse`assumptions`edgeSource`regimeApplicability`riskMemoryKey`govHypoId`owner`asOf!(
    `card_pnlAttribution;`attribution;`pnlAttribution;`v1;
    "Decompose a replay run's realized PnL into level/slope/curvature/carry/residual (R10 shock basis) + the bucketed curve risk - the quantitative 'name which edge'";
    "the level/slope/curvature axes ARE R10's parallel/slope/butterfly shock operators (not reinvented); the residual is the PLUG (realized - the four) and absorbs the unexplained curve move + costs; the carry attribution depends on R10's rollYield; reconciles to the realized total by construction";
    `na;
    "any (a PnL-decomposition capability, regime-agnostic)";
    `na;`;"qFDM evidence layer";2026.05.31);

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
