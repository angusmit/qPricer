# core/ — Foundation layer (ARCHITECTURE.md §1)

## Purpose
Math, linear algebra, stats, RNG, numerical helpers, the library loader, the config loader, and shared infrastructure (result/timing/cache/test/regression). The base of the stack.

## Dependencies
None — core depends on nothing above it. Every other layer depends on core.

## Modules
- `init.q` — the library loader: loads every layer in dependency order, sets `.qfdm.loaded`/`.qfdm.version`.
- `cfg.q` — the configuration loader: populates `.cfg` from `config/base.q` + optional `config/{env}.q` (selected by `QPRICER_ENV`). Loaded first.
- `config.q` — `.config.*`: FDM pricing-config defaults + validation (reads `.cfg.fdm`/`.cfg.iv`).
- `utilities.q` — `.utilities.*`: assertions, key/require helpers, shared numeric utilities.
- `montecarlo.q` — `.montecarlo.*`: normal/antithetic/correlated path generation, MC config + pricing-from-payoffs.
- `correlation.q` — `.correlation.*`: correlation matrices, Cholesky (`__cholesky`).
- `convergence.q` — `.convergence.*`: grid-convergence diagnostics.
- `pathdiag.q` — `.pathdiag.*`: Monte Carlo path diagnostics.
- `result.q` / `timing.q` / `cache.q` / `regression.q` / `testutil.q` — `.result`/`.timing`/`.cache`/`.regression`/`.testutil`: result wrappers, timing, memoization, regression capture, and the test-assertion helpers (`.testutil.assertTrue`/`assertNear`).

## Key API
`.cfg.*` (populated at load), `.config.defaultPricingConfig`/`defaultImpliedVolConfig`, `.montecarlo.defaultMcConfig`/`generateNormalPaths`, `.correlation.__cholesky`, `.testutil.assertTrue`/`assertNear`.

## Notes
- The config loader runs FIRST so `.cfg` is populated before any module reads it; modules read `.cfg` lazily inside their default-config functions.
- `.config` (FDM config validation) and `.cfg` (the value namespace) are distinct: `.config` returns dicts that now read from `.cfg.fdm`/`.cfg.iv`.
