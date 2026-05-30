# services/ — IPC services layer (optional, reserved, ARCHITECTURE.md §7)

## Purpose
IPC processes for always-on / multi-core scale: a gateway, an always-mapped HDB service, and a compute worker pool.

## Dependencies
Would sit beside `apps/`, talking to `data/` (HDB) and `backtest/` over IPC. Nothing depends on it.

## Modules
None — reserved skeleton.

## Key API
None yet.

## Notes
**Reserved / empty skeleton, optional and deferred.** Per ARCHITECTURE.md §7, for a single user a partitioned/splayed HDB + `peach` + `-s` slaves in one process already delivers most of the speed; the service mesh is built only when reload time or a single-process compute ceiling becomes a real pain. Over-engineering a gateway before it is needed is the failure mode.
