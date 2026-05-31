/ templates/factor_rv_replay.q - replay-mode factor relative-value runner (.template.factorRvReplay.*) (v0.82, Post-foundation R17)
/ ----------------------------------------------------------------------------
/ R17: re-run R8's factorRelativeValue strategy (curve PCA -> fade the residual) through the REALISTIC
/ end-to-end foundation R16 assembled (.workflow.runReplay: pre-register -> carded gating -> replay R12 ->
/ evidence audit R13 -> gov gates -> regime skeptic -> attribution R15 -> human packet). It REUSES R8's
/ signal (.factor.* PCA + R8's residual-fade rule) computed CAUSALLY per step, trades the ACTUAL contract
/ R11/the curve names at the residualMaturity rank, books NET PnL via .exec.fill, and emits an R12-shape
/ run record. Same SHAPE as R16's __replay so R13's audit + R15's attribution consume it.
/ -
/ THE POINT (read first): this is a VALIDATION, not a hunt for edge. R8 honestly found NO tradeable edge
/ on the OLD vectorised, frictionless, FULL-SAMPLE-PCA engine; the realistic path can only make the bar
/ HARDER (realistic costs + the CAUSAL PCA + deflation + the sealed holdout all SUBTRACT). The honest
/ expectation is that realistic evidence CONFIRMS or STRENGTHENS the no-edge finding - and that is SUCCESS.
/ -
/ DO NOT TUNE: the params are R8's, PRE-REGISTERED, read from .cfg.factor + .cfg.templates.rv (k /
/ residualMaturity / nMaturities for the PCA; lookback / entryZ / notional / txnCostRate for the fade).
/ NO new tuned config. This file does NOT edit R8's templates/factor_relative_value.q, its gates, its old
/ demo, or its tests (frozen - the historical record of the old-engine verdict); it REUSES them.
/ -
/ CAUSAL PCA - the only look-ahead in R8 was the FULL-SAMPLE PCA (it saw future curve moves). The fix:
/ at each step compute the PCA on a TRAILING lookback window of as-of curve changes (date<=asOf) ONLY -
/ NEVER the full sample. The residual-fade z (R8's rolling z, reused) is already causal. A planted FUTURE
/ curve move cannot change a past signal (the point-in-time guarantee inherited from R9's door / R11/R14).
/ -
/ ACTUAL CONTRACT: trade the contract at the residualMaturity RANK (the 3rd-nearest at rank 2) - it is
/ comfortably live (expiry >> asOf, so R13's universe check passes) and rolls FORWARD monotonically as the
/ front rolls off (so R13's rollRespected passes). Book its OWN price change (difference method, 0 on a
/ roll day - no roll gap), exactly as R12. The factor-RV SIGNAL says when/which-way; the actual contract
/ supplies the realised PnL the attribution (R15) then explains (structural slope/curvature, or residual?).
/ -
/ RESERVED-NAME NOTE: `asOf` NEVER `asof`; `comm` NOT `commodity`; builtins not shadowed (dev/avg/var/
/ first/last/signum/sums/deltas/prev/differ/inv/mmu/wsum/sublist). Right-associativity parenthesised.
/ ----------------------------------------------------------------------------

/ Merged R8 config (the PRE-REGISTERED params): the factor (PCA) params + the RV (fade/exec) params.
/ The two dicts share no keys, so the merge is clean. NO new tuned key is introduced.
.template.factorRvReplay.__cfg:{[overrides] (.cfg.factor,.cfg.templates.rv),$[0=count overrides; ()!(); overrides]};

/ The constant-maturity-by-RANK panel WITH the rank-residualMaturity contract id + price. Per date sort by
/ tenor, take the front nMaturities prices (the SAME front->deferred level panel R8's .factor.curvePanel
/ builds - NO rollBuffer, so the PCA structure is byte-faithful to R8), and ALSO pick the contractYM +
/ price at the residualMaturity rank (the ACTUAL contract we trade). Drops dates with < nMaturities live
/ contracts (as .factor.curvePanel does). PURE w.r.t. the curve history passed in.
.template.factorRvReplay.__panel:{[ch;cfg]
    nMat:cfg`nMaturities; rm:cfg`residualMaturity;
    s:`asofDate`tenor xasc ch;
    g:0!select prices:price, yms:contractYM by asofDate from s;
    g:select from g where nMat<=count each prices;
    `dates`levels`rmYM`rmPx!(g`asofDate; nMat sublist/: g`prices; (g`yms)@\:rm; (g`prices)@\:rm)
 };

/ The CAUSAL per-step residual at the residualMaturity rank: at each change-row i, run R8's PCA on the
/ TRAILING lookback window of change-rows ending at i (date<=asOf by construction), and take the residual
/ of the LAST row at the residualMaturity column. NEVER the full sample. 0 until enough rows for a k-factor
/ fit. Reuses .factor.changePanel + .factor.pca UNCHANGED. PURE + deterministic.
.template.factorRvReplay.__causalResidual:{[levels;cfg]
    chg:.factor.changePanel levels;                          / (T-1) x M curve-change panel
    nC:count chg;
    rm:cfg`residualMaturity; k:cfg`k; lb:cfg`lookback; minRows:2+k;
    {[chg;cfg;rm;k;lb;minRows;i]
        lo:0|1+i-lb;                                         / trailing lookback window (causal: rows <= i)
        win:(lo;1+i-lo) sublist chg;
        $[minRows>count win; 0f;
            (.factor.pca[win;k;cfg])[`residuals][(count win)-1; rm]]
     }[chg;cfg;rm;k;lb;minRows] each til nC
 };

/ R8's residual-fade signal, made causal. The traded series = the CUMULATIVE causal residual (R8's
/ devSrs); z = R8's causal rolling z (reuse .template.rv.__rollMean/__rollSd); position = fade extremes
/ (-sign z when |z|>entryZ), R8's exact rule. Returns the FULL-history series (sliced to a window later).
.template.factorRvReplay.__signal:{[levels;cfg]
    res:.template.factorRvReplay.__causalResidual[levels;cfg];
    devSrs:sums res;                                         / cumulative causal residual (R8's traded level)
    w:cfg`lookback; entryZ:cfg`entryZ; notional:cfg`notional;
    rm:.template.rv.__rollMean[w;devSrs]; rsd:.template.rv.__rollSd[w;devSrs];
    z:?[rsd>0f; (devSrs-rm)%rsd; 0f];
    pos:notional*?[(abs z)>entryZ; neg signum z; 0f];        / R8's fade-extremes rule
    `dev`z`pos!(devSrs;z;pos)
 };

/ The end-to-end as-of REPLAY: reuse the causal signal, trade the ACTUAL rank-rm contract, book its OWN
/ price change (difference method, 0 on a roll), cost via .exec.fill, emit an R12-shape run record
/ (`meta`steps`rollEvents`provenance) - the shape R13's audit + R15's attribution consume. Signal is
/ computed on the FULL history (causal, split-independent) then the OUTPUT is sliced to [fromDate;toDate].
.template.factorRvReplay.__replay:{[comm;fromDate;toDate;cfg]
    ch0:.data.hdb.curveHistory[comm; .data.hdb.dates comm];
    pan:.template.factorRvReplay.__panel[ch0;cfg];
    sig:.template.factorRvReplay.__signal[pan`levels;cfg];
    / change-row alignment: .factor.changePanel drops the seed, so drop the first date/contract too.
    cDates:1_ pan`dates; cYM:1_ pan`rmYM; cPx:1_ pan`rmPx;
    pos:sig`pos;                                             / one position per change-row
    allDates:cDates;
    fromD:$[null fromDate; min allDates; fromDate];
    toD:$[null toDate; max allDates; toDate];
    keep:allDates within (fromD;toD);
    dates:allDates where keep; ac:cYM where keep; px:cPx where keep; posK:pos where keep;
    n:count dates;
    / the contract's OWN price change (difference method): 0 on a roll day (contract identity changed) and
    / 0 on day 0 of the kept window (each window starts flat, the R16 convention).
    isRoll:differ ac; isRoll[0]:0b;
    dPriceFull:deltas px; dPrice:?[isRoll; 0f; dPriceFull]; dPrice[0]:0f;
    posLag:0f,-1_posK;                                       / position held INTO each day's move
    positionPnl:posLag*dPrice;
    order:deltas posK; filledQty:order;                      / position changes (full fill)
    execCfg:.exec.__resolve cfg;
    cost:{[execCfg;ord] (.exec.fill[ord; `refPrice`currentPos!(1f;0f); execCfg])`totalCost}[execCfg] each order;
    stepPnl:positionPnl-cost;
    steps:flip `stepIndex`stepDate`asOf`activeContract`isRoll`rollFrom`position`frontReturn`positionPnl`proportionalCost`slippageCost`fixedCost`rollCost`financingCost`totalCost`stepPnl`cumPnl`filledQty`refPrice`provDateTo!(
        til n; dates; dates; ac; isRoll; ?[isRoll; prev ac; 0Nj]; posK; dPrice; positionPnl;
        cost; n#0f; n#0f; n#0f; n#0f; cost; stepPnl; sums stepPnl; filledQty; px; dates);
    rollEvents:select date:stepDate, commodity:comm, fromContract:rollFrom, toContract:activeContract, reason:`rankRoll from steps where isRoll;
    `meta`steps`rollEvents`provenance!(
        `strategy`commodity`fromDate`toDate`nSteps`totalPnl!(
            `factorRelativeValue; comm; min dates; max dates; n; sum stepPnl);
        steps; rollEvents; select stepDate, asOf, provDateTo from steps)
 };

/ The runner seam for the gov gates: run the replay over [from;to] -> (date;pnl). Reuses __replay.
.template.factorRvReplay.runner:{[comm;cfg] {[comm;cfg;from;to] r:.template.factorRvReplay.__replay[comm;from;to;cfg]; select date:stepDate, pnl:stepPnl from r`steps}[comm;cfg]};

/ The replay-mode run. inputs: commodity (+ optional overrides / dateFrom / dateTo). Returns
/ `pnl`validation`meta. validation carries R8's OWN gates REUSED VERBATIM (factor-stability +
/ residual-stationarity); meta carries the R12-shape run record (the artifact .workflow.runReplay audits
/ + attributes). NOTE: this is a thin replay-mode RUNNER reusing R8's card (factorRelativeValue) for
/ gating - it is deliberately NOT registered as a new R2 template (no new card, no audit surface).
.template.factorRvReplay.run:{[inputs]
    cfg:.template.factorRvReplay.__cfg[$[`overrides in key inputs; inputs`overrides; ()!()]];
    comm:$[`commodity in key inputs; inputs`commodity; `CRUDE];
    fromD:$[`dateFrom in key inputs; inputs`dateFrom; 0Nd];
    toD:$[`dateTo in key inputs; inputs`dateTo; 0Nd];
    ch0:.data.hdb.curveHistory[comm; .data.hdb.dates comm];
    pan:.template.factorRvReplay.__panel[ch0;cfg];
    X:.factor.changePanel pan`levels;                        / the change panel R8's gates inspect
    sig:.template.factorRvReplay.__signal[pan`levels;cfg];
    stab:.template.factorRv.stability[X;cfg`k;cfg];          / R8's factor-stability gate (REUSED)
    stat:.template.rv.stationarity sig`dev;                  / R8's residual-stationarity gate (REUSED)
    gatesPass:(stab`stable) and stat`stationary;
    run:.template.factorRvReplay.__replay[comm;fromD;toD;cfg];
    pnl:select date:stepDate, pnl:stepPnl from run`steps;
    `pnl`validation`meta!(
        pnl;
        `factorStability`stationarity`shapeGate!(stab; stat;
            `passed`detail!(gatesPass; "factor-stability cos=",(string stab`cos),", residual AR1=",(string stat`ar1),", stationary=",(string stat`stationary),", gatesPass=",string gatesPass));
        `shape`commodity`nDays`ar1`halfLife`edgeSource`runRecord!(
            `factorRelativeValueReplay; comm; count pnl; stat`ar1; stat`halfLife; `structural; run))
 };

-1 "factor_rv_replay.q loaded - .template.factorRvReplay.* (R17: R8's factor-RV signal re-run causally through the realistic foundation; reuses R8's card + gates, R8 frozen)";
