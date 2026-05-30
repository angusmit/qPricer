# scripts/ — Build / CI / pipeline layer (ARCHITECTURE.md §8)

Build, CI gates, and the scheduled nightly ingest -> calibrate -> backtest -> report pipeline.

- `ingest_hdb.q` (v0.59) — builds the splayed HDB at `.cfg.paths.hdb` from the real Barchart CSVs under `.cfg.paths.barchartRoot`, reusing the parser. Idempotent: `q scripts/ingest_hdb.q`. Touches real gitignored data; not a test.
