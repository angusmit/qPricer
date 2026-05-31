/ season/season.q - curve/spread seasonality (.season.*) (v0.79, Research OS R14)
/ ----------------------------------------------------------------------------
/ ARCHITECTURE.md Part II 11.8(f) / 13-R14. An evidence-layer FEATURE capability: the curve/spread
/ SEASONALITY the first real strategy (R16) is built on - the seasonally-adjusted same-calendar-month
/ z-score IS the signal for R16's calendar-spread mean-reversion. Computed CAUSALLY through R9's door,
/ so the seasonal statistics are point-in-time by construction.
/ -
/ THE KEY POINT-IN-TIME DISCIPLINE: the same-month z-score on date D uses ONLY same-calendar-month
/ observations UP TO D (trailing/causal), obtained through R9's door (.state.asof returns date<=asOf).
/ A FULL-SAMPLE seasonal statistic is look-ahead (it standardises a past observation by future
/ same-months) - the silent killer of seasonal backtests. Below a minimum same-month N (.cfg.season)
/ the z is null + lowConfidence rather than a spurious value (thin same-month history early in the data).
/ -
/ DISTINCT from the existing signals/ seasonality (`.commodity.seasonality.*`, the alpha-library higher
/ in the stack) - .season.* is the curve/spread seasonality and does NOT touch or collide with it.
/ -
/ LOW evidence-tier: loads AFTER curve/ (reuses R10's front/deferred/slope CONVENTION + deferredIdx)
/ and reads state/ (R9) downward; opens nothing at import. ADDITIVE: reads + computes, edits nothing;
/ does NOT rewire .state.build. Registered as an R2 `season capability; carded.
/ -
/ RESERVED-NAME NOTE: parameter `asOf` NEVER `asof`; `comm` NOT `commodity`; calendar month via
/ `(`month$date) mod 12` (the `.mm` accessor is ATOM-ONLY - it errors on a date VECTOR); builtins not
/ shadowed (avg/dev/var/med/first/last/sum/group/where).
/ qSQL select/exec on a LOCAL table inside a lambda throws 'assign (skill #24) - we use mask+where+group.
/ ----------------------------------------------------------------------------

.season.defaultConfig:{[] .cfg.season};

/ The causal per-date front-deferred spread series from the as-of slice (date<=asOf), built in ONE
/ pass with plain vector ops (NO qSQL on the local slice). Per trading date: the live contracts
/ (expiry on/after the date) sorted by expiry -> front = idx0, deferred = idx deferredIdx (clamped);
/ spread = front - deferred; slope = (deferred-front)/front (R10's convention); frontYM the front's
/ contract month. Causal by construction (the slice is date<=asOf).
.season.__series:{[asOf;comm;deferredIdx]
    slice:(.state.asof[asOf;comm])`data;
    if[0=count slice; '"season: no as-of data for ",string comm];
    sl:`date`expiry xasc slice where (slice`expiry)>=slice`date;     / live rows, date then ascending expiry
    idxByDate:group sl`date;
    px:`float$sl`settle; yms:sl`contractYM;
    grp:value idxByDate;
    frontIdx:first each grp;
    defrdIdx:{[idx;k] idx[((count idx)-1)&k]}[;deferredIdx] each grp;
    frontPx:px frontIdx; defrdPx:px defrdIdx; frontYM:yms frontIdx;
    ([] date:key idxByDate; frontPx:frontPx; deferredPx:defrdPx; frontYM:frontYM;
        spread:frontPx-defrdPx; slope:?[frontPx=0f; 0n; (defrdPx-frontPx)%frontPx])
 };

/ z-score of `cur` vs the distribution `hist`, with a min-N + zero-dispersion guard (null below minN).
.season.__z:{[cur;hist;minN]
    n:count hist;
    sd:$[1<n; dev hist; 0f];
    $[(n>=minN) and sd>0f; (cur-avg hist)%sd; 0n]
 };

/ The seasonal feature dict for (asOf, comm): same-calendar-month z (the R16 signal), same-contract-
/ month z, the calendar-month seasonal factor, the deseasonalised level, the seasonally-adjusted slope,
/ + nSameMonth / lowConfidence. All CAUSAL (same-calendar-month observations up to asOf).
.season.features:{[asOf;comm]
    cfg:.cfg.season; minN:cfg`minObs;
    ser:.season.__series[asOf;comm;cfg`deferredIdx];
    effDate:max ser`date;
    cur:first ser where ser[`date]=effDate;                          / the row on the effective date
    cm:(`month$ser`date) mod 12;                                     / calendar month (0-11) of each obs (.mm is atom-only)
    curCm:(`month$effDate) mod 12;
    smMask:cm=curCm;                                                 / same calendar month, up to asOf (causal)
    smSpread:(ser`spread) where smMask;
    smSlope:(ser`slope) where smMask;
    smFront:(ser`frontPx) where smMask;
    nSameMonth:count smSpread;
    / same CONTRACT-delivery month (the front's delivery month, e.g. all March contracts historically)
    scmMask:((ser`frontYM) mod 100)=(cur`frontYM) mod 100;
    scmFront:(ser`frontPx) where scmMask;
    seasonalFactor:$[nSameMonth>=minN; avg smFront; 0n];
    `asOf`commodity`effectiveDate`spread`frontLevel`slope`sameMonthZ`sameContractMonthZ`seasonalFactor`deseasonalised`seasonalSlope`nSameMonth`lowConfidence!(
        asOf; comm; effDate; cur`spread; cur`frontPx; cur`slope;
        .season.__z[cur`spread; smSpread; minN];
        .season.__z[cur`frontPx; scmFront; minN];
        seasonalFactor;
        $[(not null seasonalFactor) and seasonalFactor>0f; (cur`frontPx)%seasonalFactor; 0n];
        $[nSameMonth>=minN; (cur`slope)-avg smSlope; 0n];
        nSameMonth; nSameMonth<minN)
 };

/ ── register as an R2 `season capability ─────────────────────────────────────

.season.contract:`version`requiredIn`requiredOut!(
    1;
    `asOf`commodity!`date`symbol;
    `asOf`commodity`sameMonthZ`seasonalFactor`deseasonalised!`date`symbol`float`float`float);
.registry.new[`season; .season.contract];
.season.register:{[name;fn;manifest] .registry.register[`season;name;fn;manifest]};
.season.get:{[name] .registry.get[`season;name]};
.season.list:{[] .registry.list `season};
.season.conforms:{[name] .registry.conforms[`season;name]};

.season.register[`curveSeasonality; .season.features;
    `contractVersion`description`in`out!(
        1;
        "causal curve/spread seasonality: same-calendar-month + same-contract-month z, seasonal factor, deseasonalised level, seasonal slope";
        `asOf`commodity!`date`symbol;
        `asOf`commodity`sameMonthZ`seasonalFactor`deseasonalised!`date`symbol`float`float`float)];

-1 "season.q loaded - .season.* causal curve/spread seasonality ready (distinct from signals/ seasonality; not opened)";
