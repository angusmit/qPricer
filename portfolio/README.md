# portfolio/ — Portfolio layer (reserved, ARCHITECTURE.md §1, §9 step 5)

## Purpose
Allocation / portfolio optimization across strategies (capital allocation, risk budgeting).

## Dependencies
Will depend on `backtest/` (strategy returns) and `analytics/` (risk). Nothing depends on it yet.

## Modules
None — reserved skeleton.

## Key API
None yet.

## Notes
**Reserved / empty skeleton.** Stood up in migration step 5 (after the execution layer); not started. `portfolio/` here is the cross-strategy allocator — distinct from `analytics/portfolio.q`, which is per-trade batch pricing + portfolio greeks.
