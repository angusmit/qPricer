# config/ — Configuration layer (ARCHITECTURE.md §2)

Holds the `.cfg` value files. The loader `core/cfg.q` (run first by `core/init.q`) populates the `.cfg` namespace:

1. `base.q` — **all** baseline defaults as domain-keyed q dicts (`.cfg.fdm`, `.cfg.iv`, `.cfg.mc`, `.cfg.calib.*`, `.cfg.analytics.*`, `.cfg.paths`). Every value equals the literal it replaced in its module, with types and key-order preserved, so loading base alone reproduces prior behavior byte-identically.
2. `{env}.q` — optional per-environment overrides, loaded **after** `base.q` only when `QPRICER_ENV` is set to that env (`dev`/`test`/`prod`). Upsert a subset; everything else falls through to base. `dev.q` is a demo override.

Unset `QPRICER_ENV` => base only => default/test behavior. An unknown env falls back to base with a warning.

Tabular reference data (calendars, instrument specs) belongs here too, as CSV/keyed tables — none yet.

**Step 2b:** the `backtest/` layer's per-strategy `defaultConfig` dicts are not yet externalized into `.cfg`.
