# analytics/ — Analytics layer (ARCHITECTURE.md §1)

## Purpose
Risk, scenarios, VaR, P&L attribution, model-risk limits, portfolio pricing, dashboards, and batch/stress/perf orchestration + reporting — everything built on top of the pricers.

## Dependencies
Depends on `core/`, `models/`, `calibration/`. Called by `apps/` and the orchestration runners.

## Modules
- Risk / scenarios: `risk.q` (`.risk` scenario reports), `riskdist.q` (`.riskdist`), `var.q` (`.var` VaR / expected shortfall), `histscen.q` (`.histscen` historical shocks), `replay.q` (`.replay` historical replay), `stress.q` (`.stress`).
- Limits: `limits.q` / `limitcheck.q` / `limitreport.q` (`.limits`/`.limitcheck`/`.limitreport`).
- Commodity model-risk: `commodityModelReport.q` (`.commodity.modelreport.*` — cross-model comparison, disagreement, model-risk limits, portfolio dashboards).
- Portfolio / P&L / reporting: `portfolio.q` (`.portfolio` batch pricing + portfolio greeks/scenarios), `pnl.q` (`.pnl`), `report.q` (`.report` CSV export), `audit.q` (`.audit`), `dashboard.q`/`dailyrisk.q` (`.dashboard`/`.dailyrisk` orchestration).
- Batch / perf: `batch.q` (`.batch`), `perfdiag.q`/`perfopt.q` (`.perfdiag`/`.perfopt`).

## Key API
`.risk.generateScenarioReport`, `.var.*`, `.portfolio.priceTrades`/`calculatePortfolioGreeks`, `.commodity.modelreport.runComparisonRisk`/`runPortfolioRiskWithLimits`, `.report.exportCsv`.

## Notes
- Model-report thresholds (greek-bump / disagreement / model-risk-limit) are externalized to `.cfg.analytics.*`.
- Per-trade isolation: batch/portfolio code wraps each trade so one bad trade returns `status=\`ERROR` rather than crashing the book.
- Output CSV directories come from caller args (`.cfg.paths.outputDir` is the documented default).
