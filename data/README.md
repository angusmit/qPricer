# data/ — Data & ingestion layer (ARCHITECTURE.md §3-4)

Barchart CSV ingestion (`parser.q`), futures-curve construction, real-data replay, and the splayed HDB load/query layer (`hdb.q`, `.data.hdb.*`, v0.59). Raw CSVs under `data/barchart/` (gitignored) are an ingestion source, not the store; the built HDB lives at `.cfg.paths.hdb` (`data/hdb`, gitignored). `hdb.q` reuses the parser's expiry/firstDate derivation and mirrors its curve-output schema, so HDB queries are byte-identical to the parser; it never opens the HDB at import.
