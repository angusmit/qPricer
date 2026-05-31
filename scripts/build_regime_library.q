\l core/init.q
/ ============================================================================
/ scripts/build_regime_library.q - build the regimeEpisodes splay alongside
/ futures / regimes in the HDB.
/ ----------------------------------------------------------------------------
/ Research OS R4 (ARCHITECTURE.md 11.3). Opens the real splayed HDB, and for each
/ curated episode in .cfg.regime.episodes computes its DOMINANT regime fingerprint
/ from the `regimes` table over the episode date range, then writes the `regimeEpisodes`
/ table splay under .cfg.paths.hdb (mirrors scripts/build_regimes.q: .Q.en + set, via
/ get/`set`, no \l). Only episodes whose ranges fall inside HDB coverage are written
/ (out-of-data lessons stay in docs/REGIME_LIBRARY.md as narrative). Idempotent.
/ Touches real, gitignored data and is NOT a test. Run AFTER scripts/build_regimes.q.
/ -
/ Usage:  q scripts/build_regime_library.q
/ ============================================================================
hdbPath:.cfg.paths`hdb;
if[0=count key hsym `$hdbPath,"/sym";
    -1 "build_regime_library: no HDB at ",hdbPath," (run scripts/ingest_hdb.q first).";
    exit 0];
if[0=count key hsym `$hdbPath,"/regimes/";
    -1 "build_regime_library: no regimes table at ",hdbPath," (run scripts/build_regimes.q first).";
    exit 0];
.regime.open hdbPath;
-1 "build_regime_library: computing dominant fingerprints for ",(string count .cfg.regime`episodes)," curated episodes ...";
t0:.z.p;
res:.regime.library.buildTable[hdbPath];
elapsed:`long$(.z.p-t0)%1000000;
-1 "build_regime_library: DONE - ",(string res`rows)," episodes -> ",hdbPath,"/regimeEpisodes/ in ",(string elapsed)," ms (",(string res`skipped)," skipped: out of HDB coverage).";
exit 0;
