\l core/init.q
/ ============================================================================
/ scripts/build_model_cards.q - persist the curated model cards (.cfg.cards) to the
/ `modelCards` splay alongside the other HDB tables.
/ ----------------------------------------------------------------------------
/ Research OS R5 (ARCHITECTURE.md knowledge plug-ins). Writes the curated card spec to
/ data/hdb/modelCards/ (mirrors scripts/build_regime_library.q: .Q.en + set, via get/`set`,
/ no \l). Idempotent. The cards are CURATED content (no computation), so this just
/ enumerates the syms and writes the splay. Touches real, gitignored data; NOT a test.
/ The knowledge plug-in also works in-memory off .cfg.cards without this step; building
/ persists it so .cards.open[hdbPath] can load it in a fresh process.
/ -
/ Usage:  q scripts/build_model_cards.q
/ ============================================================================
hdbPath:.cfg.paths`hdb;
if[0=count key hsym `$hdbPath,"/sym";
    -1 "build_model_cards: no HDB at ",hdbPath," (run scripts/ingest_hdb.q first).";
    exit 0];
t0:.z.p;
nCards:.cards.buildTable[hdbPath];
elapsed:`long$(.z.p-t0)%1000000;
-1 "build_model_cards: DONE - ",(string nCards)," cards -> ",hdbPath,"/modelCards/ in ",(string elapsed)," ms.";
exit 0;
