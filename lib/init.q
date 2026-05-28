/ init.q - silent library loader
/ Core
\l lib/utilities.q
\l lib/config.q
\l lib/product.q
\l lib/market.q
\l lib/marketbook.q
\l lib/model.q
\l lib/grid.q
\l lib/payoff.q
\l lib/boundary.q
\l lib/solver.q
\l lib/engine.q
\l lib/greeks.q
\l lib/validation.q
/ Equity models
\l lib/iv.q
\l lib/surface.q
\l lib/risk.q
\l lib/american.q
\l lib/montecarlo.q
\l lib/asian.q
\l lib/correlation.q
\l lib/basket.q
\l lib/pathdiag.q
\l lib/lookback.q
\l lib/variance.q
\l lib/convergence.q
\l lib/heston.q
\l lib/objective.q
\l lib/calibration.q
\l lib/sabr.q
\l lib/merton.q
\l lib/bates.q
/ Model quality
\l lib/limitcheck.q
\l lib/modelcheck.q
\l lib/modelcompare.q
\l lib/calibreport.q
/ Risk engine
\l lib/riskdist.q
\l lib/var.q
\l lib/histscen.q
\l lib/replay.q
\l lib/limits.q
\l lib/limitreport.q
/ Asset class and commodity
\l lib/assetclass.q
\l lib/commodityFutures.q
\l lib/commodityBlack76.q
\l lib/schwartz.q
\l lib/commoditySpread.q
\l lib/electricity.q
\l lib/commodityBacktest.q
/ Dashboard and orchestration
\l lib/dashboard.q
\l lib/dailyrisk.q
/ Portfolio and reporting
\l lib/portfolio.q
\l lib/report.q
\l lib/pnl.q
\l lib/audit.q
\l lib/regression.q
\l lib/batch.q
\l lib/result.q
\l lib/timing.q
\l lib/testutil.q
\l lib/stress.q
\l lib/perfdiag.q
\l lib/cache.q
\l lib/perfopt.q
/ Parsers and backtest
\l lib/parser.q
\l lib/backtest.q

.qfdm.loaded:1b;
.qfdm.version:"0.33";
