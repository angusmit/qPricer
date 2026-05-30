/ config/dev.q - dev environment overrides (v0.57, ARCHITECTURE.md section 2)
/ ----------------------------------------------------------------------------
/ Loaded by core/cfg.q AFTER config/base.q ONLY when QPRICER_ENV=dev. Demonstrates
/ the override mechanism: upsert a subset of base values for fast local iteration.
/ This file is NOT loaded on the default / test path, so the suite stays
/ byte-identical; it only takes effect when a developer sets QPRICER_ENV=dev.
/ ----------------------------------------------------------------------------

/ Coarser FDM grid for faster local pricing while iterating (not for production).
.cfg.fdm[`numberOfSpotSteps]:80;
.cfg.fdm[`numberOfTimeSteps]:500;

/ Smaller Monte Carlo path count for quick local runs.
.cfg.mc[`pathCount]:10000;

/ Example of repointing a data path for a local workspace layout.
.cfg.paths[`outputDir]:"output/dev";

/ Example strategy-knob override: widen the gamma-scalp delta band so local dev
/ runs rebalance less often. Only takes effect under QPRICER_ENV=dev.
.cfg.strategy.gammaScalp[`deltaBand]:0.10;
