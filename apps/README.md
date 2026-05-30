# apps/ — Applications layer (ARCHITECTURE.md §1)

## Purpose
Standalone example scripts, demos, and report generators that exercise the stack end-to-end — including on real (gitignored) data via the HDB.

## Dependencies
Depends on the whole stack below it (loads `core/init.q`). Nothing depends on apps.

## Modules
- `run_commodity_demo.q` — the commodity-stack demo runner (futures curve, Black-76, spreads, electricity).
- `examples/` — 29 standalone scenario scripts, each starting `\l core/init.q` and ending `exit 0;`. Representative: `price_european_call.q`, `calculate_greeks.q`, `price_portfolio.q`, `build_vol_surface.q` (equity core); `load_crude_curve.q`, `calibrate_crude_curve.q`, `convenience_yield_series.q`, `kalman_schwartz_smith.q`, `gas_seasonality.q`, `cross_commodity_backtest.q`, `crack_spread_real.q` (commodity, HDB-backed); `execution_gross_vs_net.q` (v0.60 gross-vs-net under slippage).

## Key API
Run any example directly: `q apps/examples/<file>.q` (each loads the library and exits). The commodity demo: `q apps/run_commodity_demo.q`.

## Notes
- The only layer that touches real data directly. Real-data examples are **data-conditional**: they prefer the splayed HDB (`.cfg.paths.hdb`) when present and fall back to the parser on CSVs, skipping gracefully if neither is present.
- Real Barchart CSVs and the built HDB are gitignored; examples are not part of `tests/run_all_tests.q` (the suite is synthetic-only).
