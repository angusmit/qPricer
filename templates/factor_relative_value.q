/ templates/factor_relative_value.q - the factorRelativeValue problem template (v0.73, R8)
/ ----------------------------------------------------------------------------
/ The FIRST problem template built on top of the new .factor.* capability: it trades the
/ mean-reversion of the curve's deviation from its factor-implied shape. Compose .factor.*
/ (PCA) -> the cumulative residual at a chosen maturity (the level the curve has drifted from
/ its k-factor reconstruction) -> the existing RV mean-reversion engine (.template.rv.*) ->
/ a per-day NET PnL series. Edge source: STRUCTURAL (factor-structure reversion).
/ -
/ TWO template-specific gates that CAN FAIL (each rejecting BEFORE the universal gov gates):
/   (a) FACTOR STABILITY - the top loading is stable across the window's two halves (cosine
/       >= minLoadingStability); if the factor structure isn't stable the "deviation" is
/       meaningless. (b) RESIDUAL STATIONARITY - the traded cumulative residual mean-reverts
/       (reuses .template.rv.stationarity's AR(1) test). New code, new (non-canonical) numbers.
/ ----------------------------------------------------------------------------

/ Factor-stability gate: |cosine| of the top PCA loading on the first vs second half of the
/ window. Both loadings are sign-fixed (front-positive), so a stable structure -> cosine ~ 1;
/ a regime break that rotates the loadings -> low cosine -> fail. PURE.
.template.factorRv.stability:{[X;k;cfg]
    nObs:count X;
    h:nObs div 2;
    if[h<3; :`cos`stable!(0n;0b)];
    a:(.factor.pca[h sublist X;k;cfg])[`loadings][;0];
    b:(.factor.pca[(neg nObs-h) sublist X;k;cfg])[`loadings][;0];
    cosSim:abs (a wsum b) % (sqrt a wsum a) * sqrt b wsum b;
    `cos`stable!(cosSim; cosSim>=cfg`minLoadingStability)
 };

/ The relativeValue-on-factors run. inputs: commodity (sym) + EITHER curveHistory (table) OR a
/ prebuilt panel (T x M change matrix) + dates; optional overrides / dateFrom / dateTo.
/ Returns `pnl`validation`meta; validation carries the factor-stability + residual-stationarity gates.
.template.factorRv.run:{[inputs]
    fcfg:.cfg.factor;
    rcfg:.cfg.templates.rv,$[`overrides in key inputs; inputs`overrides; ()!()];
    bundle:$[`panel in key inputs;
        `X`dates!(inputs`panel; inputs`dates);
        [ch:inputs`curveHistory;
         if[(`dateFrom in key inputs) and `dateTo in key inputs;
             ch:select from ch where asofDate within (inputs`dateFrom;inputs`dateTo)];
         cp:.factor.curvePanel[ch;fcfg`nMaturities];
         `X`dates!(.factor.changePanel cp`levels; 1_ cp`dates)]];
    X:bundle`X; dates:bundle`dates;
    d:.factor.pca[X;fcfg`k;fcfg];
    / the traded series: the CUMULATIVE residual at the chosen maturity (the level deviation).
    devSrs:sums (d`residuals)[;fcfg`residualMaturity];
    stab:.template.factorRv.stability[X;fcfg`k;fcfg];
    stat:.template.rv.stationarity devSrs;
    gatesPass:(stab`stable) and stat`stationary;
    detail:.template.rv.__signalPnl[devSrs;dates;rcfg];
    pnl:select date, pnl from detail;
    `pnl`validation`meta!(
        pnl;
        `factorStability`stationarity`shapeGate!(stab; stat;
            `passed`detail!(gatesPass; "factor-stability cos=",(string stab`cos),", residual AR1=",(string stat`ar1),", stationary=",(string stat`stationary),", gatesPass=",string gatesPass));
        `shape`commodity`nDays`explainedVar`edgeSource!(
            `factorRelativeValue; $[`commodity in key inputs; inputs`commodity; `]; count pnl; d`explainedVar; `structural))
 };

.template.register[`factorRelativeValue; .template.factorRv.run;
    `contractVersion`description`shape`composes`in`out!(
        1;
        "factor relative-value: fade the curve's cumulative residual from its k-factor PCA shape";
        `factorRelativeValue;
        `curvePCA;
        `commodity`curveHistory!`symbol`table;
        `pnl`validation`meta!`table`dict`dict)];

-1 "factor_relative_value.q loaded - .template.factorRv.* factorRelativeValue template registered";
