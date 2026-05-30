# scripts/ — Build / CI / pipeline layer (ARCHITECTURE.md §1, §8)

## Purpose
Build, CI gates, and the scheduled ingest → calibrate → backtest → report pipeline. Mostly reserved; the HDB ingestion script is live.

## Dependencies
Depends on the whole stack (loads `core/init.q`). Nothing depends on scripts.

## Modules
- `ingest_hdb.q` — builds the splayed HDB at `.cfg.paths.hdb` from the real Barchart CSVs under `.cfg.paths.barchartRoot`, auto-discovering commodity folders and reusing `.parser.futures` via `.data.hdb.ingest`. Idempotent (wipes + rebuilds). Run: `q scripts/ingest_hdb.q`.

## Key API
`q scripts/ingest_hdb.q` (build/rebuild the HDB).

## Notes
- `ingest_hdb.q` touches REAL, user-supplied, gitignored data and is **not** a test (the suite proves ingest↔query equivalence on synthetic data instead).
- CI gates and the nightly pipeline (ARCHITECTURE.md §8) are reserved for migration step 7.
