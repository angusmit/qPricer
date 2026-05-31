/ regime/analogue.q - regime/analogue library + risk memory (.regime.analogue.* /
/ .regime.library.*) (v0.67, Research OS R4)
/ ----------------------------------------------------------------------------
/ ARCHITECTURE.md Part II 11.3 / 13-R4. The "digested history" side of the research
/ OS: a curated, versioned library of NAMED historical regime EPISODES (grounded in
/ the HDB), each carrying its DRIVERS and its RISK MEMORY (the known strategy failure
/ modes - "what killed people here"), plus an ANALOGUE engine that, given a regime
/ state, returns the nearest historical episodes and their risk memory. This lets the
/ workflow ask "this state resembles WHAT?" and "what killed strategies in that state?"
/ from explicit, queryable records rather than an agent's recollection.
/ -
/ LAYER: regime/ stays a LOW layer - this file READS the `regimes` table (R1) + a new
/ `regimeEpisodes` table and reuses regime's own core (.regime.label); it MUST NOT
/ import gov/ or backtest/ (one-directional rule; gov ABOVE regime may consult the
/ analogue, never the reverse). Never opens the HDB at import; .regime.library.open
/ uses `get`, not `\l`. The distance/ranking functions are PURE (no HDB read) - testable.
/ -
/ HONEST DATA SCOPE: episodes are MATCHED only on periods present in the HDB (crude
/ ~2018.12 -> 2026). Out-of-data lessons (2008, 2014-16) live in REGIME_LIBRARY.md as
/ written narrative ("no data - narrative only"); they are NOT invented as episodes
/ with computed fingerprints. Only episodes whose ranges are inside HDB coverage are built.
/ ----------------------------------------------------------------------------

/ ── the analogue engine (pure: no HDB read) ──────────────────────────────────

/ The fingerprint axes shared by a regime state (.regime.label) and an episode row.
.regime.analogue.__pctAxes:`slopePct`volPct`volumePct;
.regime.analogue.__catAxes:`curveState`volState`liqState`rollPhase`seasonPhase;

/ Normalised distance between two regime fingerprints a and b (each a dict / table row
/ carrying the pct + categorical axes): a weighted Euclidean over the percentile axes
/ (each already in [0,1]) PLUS a per-mismatch categorical penalty over the discrete
/ states. Weights from .cfg.regime.analogue. Identical fingerprints -> 0. PURE.
.regime.analogue.distance:{[a;b]
    w:.cfg.regime`analogue;
    ps:.regime.analogue.__pctAxes;
    ws:w`wSlope`wVol`wVolume;
    diff:(a ps)-b ps;
    euclid:sqrt sum ws*diff*diff;
    cats:.regime.analogue.__catAxes;
    mism:sum not (a cats)=b cats;
    euclid+(w`catPenalty)*mism
 };

/ Rank an episode table by distance to a query state; return the n nearest with their
/ identity + drivers + risk-memory summary + the distance. PURE (the episodes are passed
/ in). nearest wraps this with the loaded library.
.regime.analogue.rank:{[state;episodes;n]
    if[0=count episodes; :episodes];
    dists:.regime.analogue.distance[state;] each episodes;
    r:update distance:dists from episodes;
    n sublist `distance xasc select episodeId,commodity,label,dateFrom,dateTo,driversKey,riskMemory,distance from r
 };

/ The n episodes nearest to a regime state, ranked by .regime.analogue.distance.
.regime.analogue.nearest:{[state;n] .regime.analogue.rank[state;.regime.library.episodes[];n]};

/ Convenience: the n episodes nearest to a (commodity,date)'s OWN regime fingerprint
/ (reuses .regime.label - causal). Reads the HDB (label + the loaded episodes).
.regime.analogue.forDate:{[commodity;date;n]
    .regime.analogue.nearest[.regime.label[commodity;date];n]
 };

/ ── the regime-episode library (the `regimeEpisodes` table) ──────────────────

/ Typed empty episode table (so a fresh process / the test suite works with no HDB).
.regime.library.__empty:{[]
    ([] episodeId:`symbol$(); commodity:`symbol$(); dateFrom:`date$(); dateTo:`date$();
        label:`symbol$(); driversKey:`symbol$(); riskMemory:();
        curveState:`symbol$(); volState:`symbol$(); liqState:`symbol$(); rollPhase:`symbol$();
        seasonPhase:`symbol$(); slopePct:`float$(); volPct:`float$(); volumePct:`float$())
 };

/ The loaded episode table (or empty if not opened / built). Lets .nearest work after
/ .regime.library.open, and lets a synthetic test set `regimeEpisodes directly.
.regime.library.episodes:{[] $[`regimeEpisodes in key `.; regimeEpisodes; .regime.library.__empty[]]};

/ The DOMINANT fingerprint over a set of `regimes` rows: modal discrete state on each
/ categorical axis + mean of each percentile axis. PURE (operates on supplied rows).
.regime.library.__fingerprint:{[rows]
    md:{first key desc count each group x};
    `curveState`volState`liqState`rollPhase`seasonPhase`slopePct`volPct`volumePct!(
        md rows`curveState; md rows`volState; md rows`liqState; md rows`rollPhase; md rows`seasonPhase;
        avg rows`slopePct; avg rows`volPct; avg rows`volumePct)
 };

/ Build the `regimeEpisodes` splay from the curated spec (.cfg.regime.episodes) by
/ computing each episode's dominant fingerprint from the loaded `regimes` table over its
/ date range, merging with the spec metadata (label / driversKey / riskMemory), and
/ writing the splay (mirrors the HDB build: .Q.en + set). Skips episodes with NO regimes
/ rows in range (out-of-coverage) - never fabricates a fingerprint. Idempotent.
/ Requires `regimes` loaded (.regime.open hdbPath) first. Touches the real HDB (build script).
.regime.library.buildTable:{[hdbPath]
    if[not `regimes in key `.; '"regime.library.buildTable: `regimes not loaded - call .regime.open first"];
    spec:.cfg.regime`episodes;
    rowFor:{[s]
        rr:select from regimes where commodity=s`commodity, date within (s`dateFrom;s`dateTo);
        if[0=count rr; :()];
        (`episodeId`commodity`dateFrom`dateTo`label`driversKey`riskMemory!(
            s`episodeId;s`commodity;s`dateFrom;s`dateTo;s`label;s`driversKey;s`riskMemory)),.regime.library.__fingerprint rr};
    rows:rowFor each spec;
    rows:rows where 0<count each rows;
    t:(upsert/)[.regime.library.__empty[];rows];
    t:.Q.en[hsym `$hdbPath;t];
    (hsym `$hdbPath,"/regimeEpisodes/") set t;
    `episodes`rows`skipped!(count t;count t;(count spec)-count rows)
 };

/ Load the persisted regimeEpisodes table (+ the sym domain) with `get` - NOT \l (no cwd
/ change). Errors cleanly if absent.
.regime.library.open:{[hdbPath]
    symPath:hsym `$hdbPath,"/sym";
    if[0<count key symPath; sym::get symPath];
    ep:hsym `$hdbPath,"/regimeEpisodes/";
    if[0=count key ep; '"regime.library.open: no regimeEpisodes at ",hdbPath," (run scripts/build_regime_library.q)"];
    regimeEpisodes::get ep;
    hdbPath
 };

-1 "analogue.q loaded - .regime.analogue.* / .regime.library.* ready (library not opened)";
