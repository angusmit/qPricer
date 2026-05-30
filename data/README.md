# data/ — Data & ingestion layer (ARCHITECTURE.md §1, §3-4)

## Purpose
Barchart CSV ingestion, futures-curve construction, real-data replay, and the splayed HDB load/query layer — the source of historical data for everything above.

## Dependencies
`parser.q` is intentionally standalone (no qFDM deps). `hdb.q` reuses `parser.q`. The rest depend on `core/`/`models/`. Used by `analytics/`, `backtest/`, `apps/`.

## Modules
- `parser.q` — `.parser.barchart` (AAPL options CSV loader/normaliser), `.parser.crude` (CRUDE futures-curve parser), `.parser.futures` (generalised to any `Day_<TAG>_*` product). Derives each contract's expiry=MAX / firstDate=MIN date from the data; builds `curveAt`/`curveHistory`.
- `hdb.q` — `.data.hdb.*`: ingestion (`ingest` — reuses `.parser.futures.loadAll`, `.Q.en` + `set` to a splayed `futures/`) and query (`open`/`curveAt`/`curveHistory`/`dates`), byte-identical to the parser by construction.
- `commodityFutures.q` — `.commodity.futures.*`: futures/forward curve construction helpers.
- `commodityBacktest.q` — commodity backtest helpers over parsed/curve data.
- `backtest.q` — `.parser.barchart` historical option-replay backtest functions.

## Key API
`.parser.crude.loadAll`/`curveAt`/`curveHistory`, `.parser.futures.loadAll[dir;tag]`, `.data.hdb.open`/`curveAt`/`curveHistory`/`dates`, `.data.hdb.ingest[csvDir;commodities;hdbPath]`.

## Notes
- Raw CSVs under `data/barchart/` and the built HDB at `.cfg.paths.hdb` (`data/hdb`) are **gitignored** — CSVs are an ingestion source, not the store.
- The splayed HDB is built by `scripts/ingest_hdb.q`. `core/init.q` loads the query layer but **never opens the HDB at import** (library-load independence); `.data.hdb.open` uses `get` (not `\l`) so it does not change the working directory.
- The parser is the CSV-read stage feeding ingestion and the synthetic equivalence test; it is unchanged by the HDB work.
