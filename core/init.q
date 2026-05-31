/ core/init.q - layered library loader (v0.56 migration step 1)
/ ----------------------------------------------------------------------------
/ Modules now live under the ARCHITECTURE.md layer tree (core / config / data /
/ models / calibration / analytics / signals / execution / backtest / portfolio).
/ Dependencies flow downward only. This loader preserves the EXACT proven load
/ ORDER of the former lib/init.q (only the file PATHS changed), so load behavior
/ - and therefore every test number - is byte-identical to v0.55. A later step
/ may regroup the load sequence strictly by layer once the dependency edges are
/ confirmed clean; this step changes locations and paths only, never code.
/ ----------------------------------------------------------------------------
/ Configuration layer (v0.57): populate .cfg (config/base.q + optional
/ config/{env}.q via QPRICER_ENV) FIRST, before any module that reads it.
\l core/cfg.q
/ Core (math / RNG / stats / infra)
\l core/utilities.q
\l core/config.q
/ Registry spine (Research OS R2: generic registry + per-kind contracts + conformance;
/ depends on nothing above base q; capability layers below register into it)
\l core/registry.q
/ Models (pricers + greeks) + pricing domain
\l models/product.q
\l models/market.q
\l models/marketbook.q
\l models/model.q
\l models/grid.q
\l models/payoff.q
\l models/boundary.q
\l models/solver.q
\l models/engine.q
\l models/greeks.q
\l models/validation.q
/ Calibration / equity models / core RNG (proven order preserved)
\l calibration/iv.q
\l calibration/surface.q
\l analytics/risk.q
\l models/american.q
\l core/montecarlo.q
\l models/asian.q
\l core/correlation.q
\l models/basket.q
\l core/pathdiag.q
\l models/lookback.q
\l models/variance.q
\l core/convergence.q
\l models/heston.q
\l calibration/objective.q
\l calibration/calibration.q
\l models/sabr.q
\l models/merton.q
\l models/bates.q
/ Model quality (calibration / analytics)
\l analytics/limitcheck.q
\l calibration/modelcheck.q
\l calibration/modelcompare.q
\l calibration/calibreport.q
/ Risk engine (analytics)
\l analytics/riskdist.q
\l analytics/var.q
\l analytics/histscen.q
\l analytics/replay.q
\l analytics/limits.q
\l analytics/limitreport.q
/ Asset class and commodity (models / data / calibration / signals)
\l models/assetclass.q
\l data/commodityFutures.q
\l models/commodityBlack76.q
\l models/schwartz.q
\l models/schwartz2.q
\l models/meanRevertingJump.q
\l analytics/commodityModelReport.q
\l models/commoditySpread.q
\l signals/seasonality.q
\l models/electricity.q
\l data/commodityBacktest.q
\l calibration/calibrateCurve.q
\l calibration/kalmanSchwartzSmith.q
/ Dashboard and orchestration (analytics)
\l analytics/dashboard.q
\l analytics/dailyrisk.q
/ Portfolio / reporting / infra
\l analytics/portfolio.q
\l analytics/report.q
\l analytics/pnl.q
\l analytics/audit.q
\l core/regression.q
\l analytics/batch.q
\l core/result.q
\l core/timing.q
\l core/testutil.q
\l analytics/stress.q
\l analytics/perfdiag.q
\l core/cache.q
\l analytics/perfopt.q
/ Data (parsers + real-data backtest + HDB query layer)
\l data/parser.q
\l data/hdb.q
\l data/backtest.q
/ State (Research OS R9: the as-of accessor + Market State; the single door to history;
/ wraps data/ downward; loads before the consumers; opens no HDB at import)
\l state/state.q
/ Curve (Research OS R10: the curve engine; consumes the state/ door downward; computes the
/ rich derived curve + shocks + immutable snapshots; opens nothing at import; does NOT rewire regime/)
\l curve/curve.q
/ Roll (Research OS R11: as-of-only roll discipline; depends on state/ door; independent of
/ curve/; reads downward only; opens nothing at import; does NOT rewire existing readers)
\l roll/roll.q
/ Regime (Research OS R1: market-state measurement; depends on data + signals)
\l regime/regime.q
/ Regime analogue library + risk memory (Research OS R4; depends on data + regime core)
\l regime/analogue.q
/ Execution (daily fill-and-cost simulation; used by the commodity backtest)
\l execution/execution.q
/ Backtest (strategy engine + strategies)
\l backtest/strategy.q
\l backtest/commodityStrategies.q
/ Portfolio (strategy allocation optimizer; consumes backtest returns)
\l portfolio/portfolio.q
/ Factor decomposition (Research OS R8: curve PCA capability; reads data/ downward; loads
/ before templates/ which compose it; registered as a `factor` kind; opens no HDB at import)
\l factor/factor.q
/ Templates (Research OS R6: problem-template abstraction; composes signals/backtest/
/ execution; registered as a `template` kind; independent of gov - loads before it)
\l templates/template.q
\l templates/relative_value.q
\l templates/factor_relative_value.q
/ Governance (Research OS R3: hypothesis registry + trials ledger + deflated-Sharpe
/ + gate cascade; HIGH layer - may import backtest/regime; never opens HDB at import)
\l gov/gov.q
/ Cards (Research OS R5: model-card knowledge plug-in; reads registries + gov + regime
/ library; a HIGH layer imported by nothing below it; never opens the HDB at import)
\l cards/cards.q
/ Registry population (Research OS R2): register the existing capabilities into the
/ spine (metadata only; runs after every capability layer exists). No compute change.
\l core/registry_populate.q
/ Workflow (Research OS R7: the single-process bounded research loop; composes
/ gov/cards/template/regime downward; the top functional layer; opens no HDB at import)
\l workflow/workflow.q

.qfdm.loaded:1b;
.qfdm.version:"0.76";
