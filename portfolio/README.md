# portfolio/ — Portfolio allocation layer (ARCHITECTURE.md §1, §9 step 5)

## Purpose
Allocate capital across strategies: given the per-strategy net-of-execution return series, compute weights under a choice of objective + constraints, evaluate the combined portfolio out-of-sample (causal walk-forward), and compare methods. The honest deliverable is the OOS comparison ("does optimization beat 1/N?").

## Dependencies
Depends on `core/` (math) and reuses `backtest/`'s index-based splits (`.strategy.commodityBT.__splits`), time-series perf (`__perf`), and row-table helper (`.strategy.__rowDictsToTable`). Consumes per-strategy return panels; changes nothing upstream.

## Modules
- `portfolio.q` — `.alloc.*`: covariance, allocation weights (all methods + constraints), causal walk-forward backtest, and the OOS method comparison.

## Key API
- `.alloc.covariance[returns;cfg]` — sample (optionally shrunk) covariance of an N×T train panel.
- `.alloc.weights[returns;method;cfg]` — weight vector for a method under the constraints.
- `.alloc.backtest[returns;method;cfg;splitCfg]` — walk-forward: weights per split on train, applied OOS; returns the combined OOS series' Sharpe/drawdown/turnover + realized avg weights.
- `.alloc.compare[returns;methods;cfg;splitCfg]` — OOS metrics across methods, ranked by OOS Sharpe, with equalWeight as baseline. **The deliverable.**

## Methods
`equalWeight` (1/N baseline), `inverseVol` (∝1/σ), `minVariance` (Σ⁻¹1 closed form / constrained), `riskParity` (ERC via cyclical coordinate descent — **the default**, covariance-only), `maxSharpe` (∝Σ⁻¹μ), `meanVariance` ((1/λ)Σ⁻¹μ). Return-using methods (maxSharpe/meanVariance) tend to overfit OOS — showing that is the point.

## Constraints (config-driven, `.cfg.alloc`)
`longOnly` (w≥0, default on), `fullyInvested` (Σw=1, default on) via capped-simplex bisection projection; `weightCap` (w_i≤cap); `turnoverPenalty` (κ, proximal shrink toward prev weights — the tractable rebalancing-cost proxy); optional covariance `shrinkage` toward the diagonal.

## Notes
- **Distinct from** `analytics/portfolio.q` (`.portfolio.*`, per-trade option-book pricing/greeks) and the v0.45 ensemble dashboard `.strategy.portfolio.*` (cross-PATH performance/correlation, not an allocator).
- **Causal / no look-ahead:** weights use only the split's train columns; index-based splits mean appending data adds splits without changing existing ones.
- riskParity is the default because it needs only the covariance (estimable), not expected returns (noisy). `apps/examples/portfolio_allocation.q` runs the OOS comparison on the real commodity strategies via the HDB.
