/ templates/relative_value.q - the relativeValue problem template (v0.70, Research OS R6)
/ ----------------------------------------------------------------------------
/ A genuinely DIFFERENT research shape from directional: a two-leg spread mean-reversion
/ template, generic over a spread definition (the demo instantiates a crude CALENDAR
/ spread - front vs a deferred tenor, same commodity/units, straight from the curve the
/ regime layer measures). It composes a NEW mean-reversion signal (fade extremes - the
/ opposite of momentum) and the existing fill model, and produces a per-day NET PnL series
/ that flows into the SAME universal gov gates. Edge source: STRUCTURAL (spread reversion).
/ -
/ TEMPLATE-SPECIFIC VALIDATION that CAN FAIL: a spread-stationarity gate (an OU / AR(1)
/ mean-reversion test). A non-stationary (random-walk) spread FAILS this gate and the RV
/ thesis is rejected BEFORE the universal gates even run - a gate the directional template
/ does not have. Same honesty discipline as the gov gates / R2 conformance / R5 audit.
/ This is NEW code producing NEW (non-canonical) numbers - it touches no pinned value.
/ INDEPENDENT of gov/ (no gov reference, even lazily).
/ ----------------------------------------------------------------------------

/ Causal trailing-window mean / stddev (no look-ahead), mirroring regime/__rollPct's shape.
.template.rv.__rollMean:{[w;v] {[w;v;i] lo:0|1+i-w; avg v lo+til 1+i-lo}[w;v] each til count v};
.template.rv.__rollSd:{[w;v] {[w;v;i] lo:0|1+i-w; win:v lo+til 1+i-lo; $[1<count win; dev win; 0f]}[w;v] each til count v};

/ Spread-stationarity (the template-specific gate). Regress spread_t on spread_{t-1}
/ (AR(1)); a mean-reverting spread has 0<=b<1 (and an OU half-life -log2/log b); a random
/ walk has b ~ 1. Stationary iff b < maxAr1 (and b > -1). PURE.
.template.rv.stationarity:{[spread]
    n:count spread;
    if[n<3; :`ar1`halfLife`stationary!(0n;0n;0b)];
    x:-1_spread; y:1_spread;
    mx:avg x; my:avg y;
    denom:sum (x-mx)*x-mx;
    b:$[denom>0f; (sum (x-mx)*y-my)%denom; 0n];
    hl:$[(not null b) and (b>0f) and b<1f; (neg log 2f)%log b; 0w];
    stationary:(not null b) and (b<.cfg.templates.rv`maxAr1) and b> -1f;
    `ar1`halfLife`stationary!(b;hl;stationary)
 };

/ The mean-reversion signal + NET PnL on a spread series. Fade extremes: position
/ = -sign(z) when |z|>entryZ (z = causal rolling z-score), held one day into the next
/ spread move; cost on each day's position CHANGE routes through the existing fill model
/ (.exec.fill, refPrice=1 / notional units, so proportional cost = rate*|dpos| and slippage
/ = |dpos|*bps). PURE w.r.t. inputs (no HDB). Returns a per-day detail table.
.template.rv.__signalPnl:{[spread;dates;cfg]
    w:cfg`lookback; entryZ:cfg`entryZ; notional:cfg`notional;
    rm:.template.rv.__rollMean[w;spread];
    rsd:.template.rv.__rollSd[w;spread];
    z:?[rsd>0f; (spread-rm)%rsd; 0f];
    pos:notional*?[(abs z)>entryZ; neg signum z; 0f];
    posLag:prev pos; posLag[0]:0f;
    dSpread:deltas spread; dSpread[0]:0f;
    gross:posLag*dSpread;
    execCfg:.exec.__resolve cfg;
    dpos:deltas pos; dpos[0]:pos 0;
    cost:{[execCfg;ord] (.exec.fill[ord; `refPrice`currentPos!(1f;0f); execCfg])`totalCost}[execCfg] each dpos;
    ([] date:dates; spread; z; pos; gross; cost; pnl:gross-cost)
 };

/ Build a calendar-spread series (front minus the deferredIdx-th tenor) from a curve
/ history (asofDate/tenor/price/...), mirroring the regime panel's front/deferred pick.
.template.rv.__spreadFromCurveHist:{[ch;deferredIdx]
    g:0!select front:first price, deferred:price@(deferredIdx)&-1+count price by asofDate from `asofDate`tenor xasc ch;
    `date`spread!(g`asofDate; (g`front)-g`deferred)
 };

/ The relativeValue template run. inputs: commodity (sym), curveHistory (table); optional
/ deferredIdx, overrides (a cfg patch), dateFrom/dateTo (the gov runner seam).
/ Returns `pnl`validation`meta (validation carries the spread-stationarity gate result).
.template.rv.run:{[inputs]
    rcfg:.cfg.templates.rv,$[`overrides in key inputs; inputs`overrides; ()!()];
    deferredIdx:$[`deferredIdx in key inputs; inputs`deferredIdx; rcfg`deferredIdx];
    ch:inputs`curveHistory;
    if[(`dateFrom in key inputs) and `dateTo in key inputs;
        ch:select from ch where asofDate within (inputs`dateFrom;inputs`dateTo)];
    sp:.template.rv.__spreadFromCurveHist[ch;deferredIdx];
    spread:sp`spread; dates:sp`date;
    stat:.template.rv.stationarity spread;
    detail:.template.rv.__signalPnl[spread;dates;rcfg];
    pnl:select date, pnl from detail;
    `pnl`validation`meta!(
        pnl;
        `stationarity`shapeGate!(stat;
            `passed`detail!(stat`stationary; "spread-stationarity gate: AR1=",(string stat`ar1),", half-life=",(string stat`halfLife),", stationary=",string stat`stationary));
        `shape`commodity`nDays`ar1`halfLife`edgeSource!(
            `relativeValue; $[`commodity in key inputs; inputs`commodity; `]; count pnl; stat`ar1; stat`halfLife; `structural))
 };

.template.register[`relativeValue; .template.rv.run;
    `contractVersion`description`shape`composes`in`out!(
        1;
        "relative-value spread mean-reversion (a new shape) with its own spread-stationarity gate";
        `relativeValue;
        `meanReversionSpread;
        `commodity`curveHistory!`symbol`table;
        `pnl`validation`meta!`table`dict`dict)];

-1 "relative_value.q loaded - .template.rv.* relativeValue template registered";
