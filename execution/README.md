# execution/ — Execution-simulation layer (ARCHITECTURE.md §1, §5)

## Purpose
Turns a strategy's target position change (an order) into a simulated fill and explicit cost components, so a backtest can book realized (net-of-execution) PnL and report GROSS vs NET performance. A daily fill-and-cost model matched to daily-settle data — not a tick / limit-order-book simulator.

## Modules
- `execution.q` — `.exec.*`: the daily fill-and-cost model (config + pure fill function).

## Dependencies
Reads `.cfg.exec` (config layer). Depends on nothing above it; the commodity backtest (`backtest/`) calls into it. Pure/stateless, so the equity engine can adopt it later (step 4b).

## Key API
- `.exec.defaultConfig[]` → the `.cfg.exec` dict.
- `.exec.__resolve[stratCfg]` → resolved exec config for a run: `.cfg.exec` + an optional per-run `exec` sub-config, with `proportionalRate` forced to the strategy's own `txnCostRate` (so the default reproduces the legacy cost).
- `.exec.fill[order;ctx;cfg]` → dict `filledQty`/`fillPrice`/`proportionalCost`/`slippageCost`/`fixedCost`/`totalCost`. `order` is the signed target-minus-current position change; `ctx` carries `refPrice`/`barVolume`/`volatility`/`currentPos`.

## Cost model
- Filled quantity: full unless `participationCap` is set and a bar volume is present → `filled = sign(order)·min(|order|, cap·barVolume)` (partial fills make the held position lag target).
- Proportional cost = `proportionalRate · |filled|` — the legacy `txnCostRate · turnover` (no price multiply; positions are vol-targeted notional scales).
- Slippage cost = `|filled| · refPrice · slippageBps/1e4` (adverse; optional vol-scaling and size/ADV impact).
- Fixed-per-trade cost when anything trades. `totalCost` = sum of the three.

## Notes
- **Frictionless by default, byte-identical.** With `.cfg.exec` defaults (zero slippage/fixed, no cap, proportional at the strategy's txnCostRate) the order fills fully and `totalCost` equals the legacy cost, so the wired backtests reproduce the pre-execution-layer numbers to the digit. Realism (slippage / fixed / participation cap / impact) is opt-in via a per-run `exec` sub-config in the strategy config.
- **Both desks fully wired (Step 4 complete).** Commodity BT (v0.60): order = vol-targeted Δposition, refPrice = front price. Equity hedge (v0.61): the shared `.strategy.__hedgeStep`/`__hedgeInit` helper passes the hedge order as dollar NOTIONAL with refPrice=1, reproducing the price-scaled legacy hedge cost `|hedgeTrade|*spot*rate`. Equity OPTION LEGS (v0.63, step 4c): every discrete leg cost (entry/roll/book/vega/lifecycle, plus spread/futures-leg and inline-hedge costs) routes through `.strategy.__legCost[notional;stratCfg]` (order = dollar notional, refPrice=1), so slippage = `|notional|*bps` is the premium-scaled option bid-ask. Cost threads into both stepPnl and cash so `deltaPV==stepPnl` holds. A per-run `exec` sub-config reaches the hedge via `.strategy.__withExec` (all hedged strategies) and the option legs directly via `__legCost`.
- The participation cap applies to the continuous hedge leg; discrete option legs fill fully (slippage/commission as cost).
- Because positions trade against price/premium-scaled slippage, even a few bps is a large drag — `apps/examples/execution_gross_vs_net.q` (commodity, HDB) and `execution_gross_vs_net_equity.q` (gamma-scalp hedge slippage + iron-condor option-leg slippage, synthetic) show the gross-vs-net decay.
