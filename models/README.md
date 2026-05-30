# models/ — Pricing layer (ARCHITECTURE.md §1)

## Purpose
Valuation + greeks for every product/model: the equity BS/FDM core, the commodity models, the stochastic-vol / jump models, the exotics, and the pricing domain (trade/market/model dicts).

## Dependencies
Depends on `core/` (and `config/` via `.cfg`). Called by `calibration/`, `analytics/`, `backtest/`.

## Modules
- Pricing domain: `product.q` (`.product` trade schema + knock-in helpers), `market.q`/`marketbook.q` (`.market`/`.marketbook` single/multi-symbol data), `model.q` (`.model` model selectors).
- FDM core: `grid.q`/`payoff.q`/`boundary.q` (`.grid`/`.payoff`/`.boundary` solver building blocks), `solver.q` (`.solver` explicit + Crank-Nicolson), `engine.q` (`.engine` — the public pricing API), `greeks.q` (`.greeks`), `validation.q` (`.validation` analytic checks), `american.q` (`.american` early-exercise boundary).
- Monte Carlo / exotics: `asian.q`, `basket.q`, `lookback.q`, `variance.q` (Asian/basket/lookback/variance pricers, MC).
- Stochastic vol / jumps: `heston.q`, `sabr.q`, `merton.q`, `bates.q`.
- Commodity: `assetclass.q` (`.assetclass` routing registry), `commodityBlack76.q` (`.commodity.black76`), `schwartz.q` (`.commodity.schwartz`, one-factor), `schwartz2.q` (`.commodity.schwartz2`, two-factor), `meanRevertingJump.q` (`.commodity.mrjump`), `commoditySpread.q` (`.commodity.spread`, Kirk/Margrabe), `electricity.q` (`.electricity`).

## Key API
`.engine.priceOption[trade;marketData;model;config]` / `.engine.priceOptionWithGrid`, `.greeks.calculateGreeks`, `.commodity.schwartz2.futuresCurve`, `.commodity.spread.kirkPrice`/`margrabePrice`, `.assetclass.__productMap`/`__modelMap`.

## Notes
- Models are caller-supplied-parameter (no library default-param dicts), so no `.cfg` sweep was needed here; MC/FDM defaults come from `.cfg.mc`/`.cfg.fdm` in core.
- Knock-in barriers are priced as `vanilla − knock-out` inside the engine; Crank-Nicolson is European-vanilla only; local vol is European-vanilla explicit-FDM only.
- When adding a product/model, update `assetclass.q`'s routing maps.
