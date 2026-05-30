/ portfolio/portfolio.q - strategy allocation optimizer (.alloc.*) (v0.62)
/ ----------------------------------------------------------------------------
/ Migration step 5 (ARCHITECTURE.md sections 6 / 9). Given the per-strategy
/ NET-of-execution return series (produced on both desks via .exec.fill), compute
/ allocation weights under a choice of objective + constraints, evaluate the
/ combined portfolio OUT-OF-SAMPLE (causal walk-forward), and COMPARE methods.
/ The honest deliverable is the OOS comparison ("does optimization beat 1/N?"),
/ not a single optimal portfolio.
/ -
/ ADDITIVE: this layer CONSUMES existing returns; it changes nothing upstream.
/ Distinct namespace `.alloc.*` - NOT analytics' option-portfolio `.portfolio.*`,
/ NOT the v0.45 cross-PATH ensemble dashboard `.strategy.portfolio.*`.
/ -
/ Convention: a return panel is a list of N strategy return series (N rows, each
/ length T) - returns[i] is strategy i's daily series. Weights are an N-vector.
/ -
/ DEFAULT METHOD = riskParity (equal risk contribution): uses only the covariance
/ (estimable), not expected returns (the noisy thing) - stable, doesn't overfit.
/ equalWeight (1/N) is always the baseline. Return-using methods (maxSharpe,
/ meanVariance) are available but expected to overfit OOS - showing that is the point.
/ -
/ CAUSALITY: weights are estimated on TRAIN columns only and applied to the OOS
/ test window; walk-forward splits are index-based over the return dates (reusing
/ .strategy.commodityBT.__splits), so appending data adds splits without changing
/ existing ones. No future return ever enters a weight.
/ ----------------------------------------------------------------------------

.alloc.defaultConfig:{[] .cfg.alloc};

/ Diagonal matrix from a vector (v on the diagonal, zeros off).
.alloc.__diagM:{[v] v*(til count v)=\:til count v};

/ Diagonal (vector) of a square matrix.
.alloc.__diag:{[m] m ./: (til count m),'til count m};

/ Sample covariance (T-1 normalised) from an N x T train panel, with optional
/ shrinkage toward the diagonal target: (1-d)*S + d*diag(S). d from cfg.shrinkage.
.alloc.covariance:{[returns;cfg]
    n:count returns;
    t:count first returns;
    mu:avg each returns;
    xc:returns-mu;
    covm:(xc mmu flip xc)%t-1;
    delta:$[`shrinkage in key cfg; cfg`shrinkage; 0f];
    if[delta>0f; covm:((1f-delta)*covm)+delta*.alloc.__diagM .alloc.__diag covm];
    covm
 };

/ Euclidean projection of v onto {w: 0<=w<=cap, sum w = 1} by bisection on the
/ clip threshold tau (sum of clipped values is monotone decreasing in tau).
.alloc.__projectCappedSimplex:{[v;cap]
    step:{[v;cap;lohi] mid:0.5*sum lohi; $[1f<sum cap&0f|v-mid; (mid;lohi 1); (lohi 0;mid)]};
    lohi:120 step[v;cap]/((min v)-cap;max v);
    tau:0.5*sum lohi;
    cap&0f|v-tau
 };

/ Risk-parity (ERC) weights via cyclical coordinate descent: minimise
/ (1/2)w'Sigma w - b*sum log w_i (b=1), whose stationarity w_i(Sigma w)_i = b gives
/ EQUAL risk contributions; normalise to sum 1. Long-only by construction.
.alloc.__riskParity:{[covm;cfg]
    n:count covm;
    dv:.alloc.__diag covm;
    w:reciprocal sqrt dv; w:w%sum w;
    maxIter:$[`rpMaxIter in key cfg; cfg`rpMaxIter; 300];
    i:0;
    while[i<maxIter;
        wOld:w;
        j:0;
        while[j<n;
            sw:covm mmu w;
            a:covm[j;j];
            bj:sw[j]-w[j]*a;
            w[j]:(neg[bj]+sqrt (bj*bj)+4f*a)%2f*a;
            j+:1];
        if[1e-12>sqrt sum (w-wOld)*w-wOld; i:maxIter];
        i+:1];
    w%sum w
 };

/ Apply weight constraints. longOnly -> project onto the capped simplex (0<=w<=cap,
/ sum 1). Else fullyInvested -> normalise to sum 1; else return the raw weights
/ (the unconstrained meanVariance case). turnoverPenalty is applied separately.
.alloc.__applyConstraints:{[raw;cfg]
    $[cfg`longOnly; .alloc.__projectCappedSimplex[raw;cfg`weightCap];
      cfg`fullyInvested; raw%sum raw;
      raw]
 };

/ Proximal turnover penalty: shrink w toward prevWeights by kappa (the tractable
/ proxy for rebalancing cost). ||blend-prev||1 = ||w-prev||1 / (1+kappa) < ||w-prev||1.
.alloc.__applyTurnover:{[w;prevW;cfg]
    kappa:$[`turnoverPenalty in key cfg; cfg`turnoverPenalty; 0f];
    if[kappa<=0f; :w];
    blended:(w+kappa*prevW)%1f+kappa;
    blended%sum blended
 };

/ Allocation weights for one method under the configured constraints.
.alloc.weights:{[returns;method;cfg]
    n:count returns;
    if[method=`equalWeight; :n#1f%n];
    covm:.alloc.covariance[returns;cfg];
    mu:avg each returns;
    sig:sqrt .alloc.__diag covm;
    lambda:$[`riskAversion in key cfg; cfg`riskAversion; 1f];
    raw:$[method=`inverseVol;   reciprocal sig;
          method=`minVariance;  (inv covm) mmu n#1f;
          method=`maxSharpe;    (inv covm) mmu mu;
          method=`meanVariance; (1f%lambda)*(inv covm) mmu mu;
          method=`riskParity;   .alloc.__riskParity[covm;cfg];
          '"alloc.weights: unknown method ",string method];
    w:.alloc.__applyConstraints[raw;cfg];
    $[`prevWeights in key cfg; .alloc.__applyTurnover[w;cfg`prevWeights;cfg]; w]
 };

/ Walk-forward OOS backtest of one allocation method: estimate weights per split on
/ that split's TRAIN columns only, apply to the OOS test columns, accumulate the
/ combined portfolio return series, and score it. Reuses the commodity backtest's
/ index-based splits (causal) and time-series perf helper.
.alloc.backtest:{[returns;method;cfg;splitCfg]
    n:count returns; t:count first returns;
    sc:.strategy.commodityBT.defaultSplitCfg[]; if[count splitCfg; sc:sc,splitCfg];
    splits:.strategy.commodityBT.__splits[sc`scheme;t;sc`trainSpan;sc`testSpan;sc`maxSplits];
    if[0=count splits; '"alloc.backtest: no valid splits for the given spans/window"];
    annDays:$[`annualizationDays in key cfg; cfg`annualizationDays; 252f];
    oosRet:`float$();
    weightsList:();
    prevW:n#1f%n;
    i:0;
    while[i<count splits;
        sp:splits i;
        trainCols:(sp`trainStartIdx)+til 1+(sp`trainEndIdx)-sp`trainStartIdx;
        testCols:(1+sp`trainEndIdx)+til (sp`testEndIdx)-sp`trainEndIdx;
        w:.alloc.weights[returns[;trainCols];method;cfg,(enlist `prevWeights)!enlist prevW];
        oosRet,:w mmu returns[;testCols];
        weightsList,:enlist w;
        prevW:w;
        i+:1];
    perf:.strategy.commodityBT.__perf[oosRet;annDays;1f];
    wMat:weightsList;
    avgW:avg wMat;
    turn:$[1<count wMat; avg sum each abs (1_wMat)-(-1_wMat); 0f];
    `method`oosSharpe`oosAnnReturn`oosAnnVol`oosMaxDrawdown`oosHitRate`avgTurnover`nSplits`oosLen`avgWeights!(
        method;perf`sharpe;perf`annualReturn;perf`annualVol;perf`maxDrawdown;perf`hitRate;turn;count splits;count oosRet;avgW)
 };

/ THE DELIVERABLE: OOS portfolio metrics across methods, ranked by OOS Sharpe, with
/ equalWeight as the baseline. Returns a table (one row per method) + realized avg weights.
.alloc.compare:{[returns;methods;cfg;splitCfg]
    rows:.alloc.backtest[returns;;cfg;splitCfg]'[methods];
    `oosSharpe xdesc .strategy.__rowDictsToTable rows
 };

-1 "portfolio.q loaded - .alloc.* strategy allocator ready (riskParity default)";
