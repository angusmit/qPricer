/ regime/regime.q - Market State Engine v1 (.regime.*) (v0.64, Research OS R1)
/ ----------------------------------------------------------------------------
/ ARCHITECTURE.md Part II 11.1 / 11.2 / 13-R1. Deterministic MEASUREMENT of a
/ per-(commodity,date) regime fingerprint from the existing HDB futures data:
/ curve state, realized-vol state, liquidity state, roll proximity, seasonal phase.
/ It does NOT predict or opine - interpretation lives in the agent layer. This is
/ the seam the rest of the research OS hangs off (regime-conditional gates, the
/ analogue engine, the regime library).
/ -
/ ADDITIVE: depends only on data/ (the HDB query layer) and signals/ (seasonality);
/ never imports backtest/ or above. .regime.breakdown (Part 3) is a COMPOSITION that
/ takes daily PnL + regime labels as plain inputs - no upward dependency. Loads at
/ import but never opens the HDB (a fresh checkout / the test suite has no HDB).
/ -
/ Pure core: __panelFromLong (build a per-date front/deferred panel from the futures
/ long table) + __labelPanel (label a panel) are HDB-free and known-answer testable.
/ .regime.series / .regime.label wrap the HDB pull; thresholds live in .cfg.regime.
/ ----------------------------------------------------------------------------

.regime.defaultConfig:{[] .cfg.regime};

/ Rolling causal percentile: for each i, the fraction of the trailing window
/ (the last `w` observations including i) that are <= v[i]. In [0,1]; uses only
/ current + past values (no look-ahead). A fraction, not a nearest-rank index.
.regime.__rollPct:{[w;v]
    f:{[w;v;i] lo:0|1+i-w; win:v lo+til 1+i-lo; (sum win<=v i)%count win};
    f[w;v] each til count v
 };

/ Build a per-date panel from a futures long table (the .data.hdb.__longTable shape:
/ contractYM,expiry,firstDate,date,settle,open,high,low,volume). Per date: the FRONT
/ contract (smallest tenor among contracts alive that day) gives frontSettle / frontVolume
/ / frontExpiry / frontTenorDays; the DEFERRED contract is the deferredIdx-th by tenor
/ (clamped to the last available). Returns one row per date, sorted by date.
.regime.__panelFromLong:{[lt;cfg]
    deferredIdx:cfg`deferredIdx;
    alive:`date`tenorDays xasc update tenorDays:`int$expiry-date from select from lt where expiry>=date;
    0!select frontSettle:first settle, frontVolume:`float$first volume, frontExpiry:first expiry,
        frontTenorDays:first tenorDays,
        deferredSettle:settle@(deferredIdx)&-1+count settle
        by date from alive
 };

/ Label a per-date panel (the __panelFromLong output) into the regime table:
/ date + the five label axes + the three underlying percentiles. Pure + deterministic.
.regime.__labelPanel:{[panel;cfg]
    dts:panel`date;
    frontSettle:panel`frontSettle;
    / curve slope (relative front-vs-deferred): <0 backwardation, >0 contango.
    slope:(`float$(panel`deferredSettle)-frontSettle)%frontSettle;
    flatThr:cfg`flatSlopeThreshold;
    curveState:?[slope<neg flatThr;`backwardation;?[slope>flatThr;`contango;`flat]];
    / realized vol from front-settle log returns, zeroed on roll days (front expiry change).
    logret:0f^deltas log frontSettle;
    rollDay:(prev panel`frontExpiry)<>panel`frontExpiry;
    ret:?[rollDay;0f;logret]; ret[0]:0f;
    rvol:(sqrt cfg`annualizationDays)*cfg[`volLookback] mdev ret;
    / causal rolling percentiles.
    pl:cfg`pctLookback;
    slopePct:.regime.__rollPct[pl;slope];
    volPct:.regime.__rollPct[pl;rvol];
    volumePct:.regime.__rollPct[pl;panel`frontVolume];
    volState:?[volPct<cfg`lowPctThreshold;`low;?[volPct>cfg`highPctThreshold;`high;`normal]];
    liqState:?[volumePct<cfg`thinPctThreshold;`thin;?[volumePct>cfg`deepPctThreshold;`deep;`normal]];
    rollPhase:?[(panel`frontTenorDays)<=cfg`rollNearDays;`near;`mid];
    / seasonal phase: REUSE signals/seasonality. fracYear encodes the calendar month.
    monthIdx0:(`mm$dts) mod 12;
    timeYears:(`year$dts)+(monthIdx0+0.5)%12f;
    seasFactor:.commodity.seasonality.factor[timeYears;cfg`seasonCfg];
    seasonPhase:?[seasFactor>cfg`seasonPosThreshold;`peak;?[seasFactor<cfg`seasonNegThreshold;`trough;`neutral]];
    ([] date:dts; curveState; volState; liqState; rollPhase; seasonPhase; slopePct; volPct; volumePct)
 };

/ Regime table over a date vector for one commodity (pulls the futures long table).
.regime.series:{[commodity;dates]
    cfg:.regime.defaultConfig[];
    lt:select from .data.hdb.__longTable[commodity] where date in dates;
    panel:.regime.__panelFromLong[lt;cfg];
    .regime.__labelPanel[panel;cfg]
 };

/ Regime labels for ONE (commodity,date) - causal: computed over that commodity's
/ history up to and including `date`, returning the row for `date`.
.regime.label:{[commodity;date]
    allDates:.data.hdb.dates commodity;
    s:.regime.series[commodity;allDates where allDates<=date];
    last select from s where date=date
 };

/ Build the regimes splay alongside futures under hdbPath (mirrors .data.hdb.ingest:
/ .Q.en + p# on commodity + set). Idempotent. Touches the real HDB (build script only).
.regime.buildTable:{[commodities;hdbPath]
    commodities:(),commodities;
    hp:hsym `$hdbPath;
    t:raze {[c] update commodity:c from .regime.series[c;.data.hdb.dates c]} each commodities;
    t:`commodity`date xasc t;
    t:.Q.en[hp;t];
    t:update `p#commodity from t;
    (hsym `$hdbPath,"/regimes/") set t;
    `rows`commodities!(count t;commodities)
 };

/ Load the persisted regimes table (and the sym domain) with `get` - NOT \l, so it
/ never changes the process working directory. Errors cleanly if absent.
.regime.open:{[hdbPath]
    symPath:hsym `$hdbPath,"/sym";
    regPath:hsym `$hdbPath,"/regimes/";
    if[0=count key regPath; '"regime.open: no regimes table at ",hdbPath," (run scripts/build_regimes.q)"];
    sym::get symPath;
    regimes::get regPath;
    hdbPath
 };

/ Part 3: regime-conditional performance breakdown. Given a (date,pnl) table and a
/ per-day regime-label table, group the PnL by the chosen axis bucket and score each
/ bucket by REUSING .strategy.commodityBT.__perf, plus a `(blended)` all-data row. Pure
/ composition - does not touch the backtest engine. The blended row reproduces the
/ backtest's headline number exactly.
.regime.breakdown:{[pnlByDate;regimeLabels;axis]
    annDays:.cfg.regime`annualizationDays;
    labTbl:`date xkey ([] date:regimeLabels`date; bucket:regimeLabels axis);
    m:select from (pnlByDate lj labTbl) where not null bucket;
    perfRow:{[annDays;sub;bname]
        p:.strategy.commodityBT.__perf[sub`pnl;annDays;1f];
        `bucket`nDays`totalPnl`annualReturn`annualVol`sharpe`maxDrawdown`hitRate!(
            bname;count sub;p`totalPnl;p`annualReturn;p`annualVol;p`sharpe;p`maxDrawdown;p`hitRate)};
    buckets:asc distinct m`bucket;
    rows:{[annDays;m;perfRow;b] perfRow[annDays;select from m where bucket=b;b]}[annDays;m;perfRow] each buckets;
    rows,enlist perfRow[annDays;m;`$"(blended)"]
 };

-1 "regime.q loaded - .regime.* Market State Engine ready (not opened)";
