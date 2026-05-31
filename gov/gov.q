/ gov/gov.q - research governance v1 (.gov.*) (v0.65, Research OS R3)
/ ----------------------------------------------------------------------------
/ ARCHITECTURE.md Part II 11.4 / 13-R3. The anti-overfitting heart: a hypothesis
/ REGISTRY (pre-registration), an append-only TRIALS LEDGER (the honest N - the
/ multiple-testing denominator), the deflated-Sharpe statistic (Bailey & Lopez de
/ Prado), and an ordered GATE CASCADE (thesis -> cost -> deflated Sharpe ->
/ walk-forward). R1 (regime/) SURFACES regime-conditional structure; R3 JUDGES it
/ and refuses to be fooled by a tempting small-sample slice.
/ -
/ HIGH layer: deps flow downward, so gov/ MAY import backtest/ (it sits ABOVE it).
/ It reuses .strategy.commodityBT.__perf (per-period <-> annualised Sharpe) and
/ __splits (causal walk-forward index splits, the same pattern .alloc.compare uses),
/ reuses .validation.__normalCdf for Phi (the BS/pricing core normal CDF, accuracy
/ 7.5e-8 - NOT reimplemented), and judges the regime-conditional results from regime/.
/ -
/ NON-OPTIONAL logging without touching engine code: .gov.run is a WRAPPER that
/ ALWAYS logs a trial per regime bucket THEN evaluates the gates - the backtest
/ ENGINE is never edited, so the full suite stays byte-identical. Logging is a
/ convention (all research goes through .gov.run), not an engine hook.
/ -
/ Pure functions (.gov.psr / .gov.deflatedSharpe / .gov.phiInv / .gov.evaluate) are
/ trivially testable on synthetic inputs; thresholds + the realistic cost config
/ live in .cfg.gov. The loader NEVER opens an HDB at import (.gov.open uses `get`,
/ not \l - same trap as .data.hdb.open / .regime.open); a fresh process gets EMPTY
/ in-memory tables via .gov.init.
/ ----------------------------------------------------------------------------

.gov.defaultConfig:{[] .cfg.gov};

/ ── the statistical core: PSR / deflated Sharpe (pure functions) ──────────────

/ Phi: standard normal CDF - REUSE the BS/pricing core (Abramowitz & Stegun, 7.5e-8).
.gov.__phi:{[x] .validation.__normalCdf x};

/ Horner evaluation of a polynomial with coefficients `co` (highest degree first)
/ at scalar z: ((..(co0*z+co1)*z+co2)..)+coN.
.gov.__horner:{[co;z] {[z;acc;c] (acc*z)+c}[z]/[0f;co]};

/ PhiInv: inverse standard normal CDF via Acklam's rational approximation
/ (relative error < 1.15e-9). Scalar; returns +/-0w at the open-interval ends.
/ Added here because no inverse-normal exists elsewhere in the codebase (the IV /
/ calibration solver root-finds the FORWARD CDF; it has no closed inverse).
.gov.phiInv:{[p]
    if[p<=0f; :neg 0w];
    if[p>=1f; :0w];
    a:-39.69683028665376 220.9460984245205 -275.9285104469687 138.3577518672690 -30.66479806614716 2.506628277459239;
    b:-54.47609879822406 161.5858368580409 -155.6989798598866 66.80131188771972 -13.28068155288572 1.0;
    c:-0.007784894002430293 -0.3223964580411365 -2.400758277161838 -2.549732539343734 4.374664141464968 2.938163982698783;
    d:0.007784695709041462 0.3224671290700398 2.445134137142996 3.754408661907416 1.0;
    plow:0.02425; phigh:1f-plow;
    $[p<plow;
        [qq:sqrt neg 2f*log p; (.gov.__horner[c;qq])%.gov.__horner[d;qq]];
      p>phigh;
        [qq:sqrt neg 2f*log 1f-p; neg (.gov.__horner[c;qq])%.gov.__horner[d;qq]];
        [qq:p-0.5; r:qq*qq; qq*(.gov.__horner[a;r])%.gov.__horner[b;r]]]
 };

/ Probabilistic Sharpe Ratio: P(true per-period Sharpe > srStar | observed srHat,
/ n returns, skewness skew, kurtosis kurt [non-excess; normal=3]).
/ PSR = Phi( (srHat-srStar)*sqrt(n-1) / sqrt(1 - skew*srHat + ((kurt-1)/4)*srHat^2) ).
.gov.psr:{[srHat;srStar;n;skew;kurt]
    if[n<=1; :0n];
    num:(srHat-srStar)*sqrt n-1f;
    t1:skew*srHat;
    t2:((kurt-1f)%4f)*srHat*srHat;
    denArg:(1f-t1)+t2;
    if[denArg<=0f; :0n];
    .gov.__phi num%sqrt denArg
 };

/ Expected MAX per-period Sharpe under the null across nTrials independent trials
/ (Bailey & Lopez de Prado): srExpMax = sqrt(V) * ((1-g)*PhiInv(1-1/N) + g*PhiInv(1-1/(N*e))).
/ V = variance of the per-period Sharpe estimates across the family's trials, g =
/ Euler-Mascheroni. N<=1 (single trial) or V<=0 -> 0 (no multiple-testing penalty).
.gov.expectedMaxSR:{[nTrials;varSR]
    if[nTrials<=1; :0f];
    if[varSR<=0f; :0f];
    g:0.5772156649015329;
    e:exp 1f;
    z1:.gov.phiInv 1f-1f%nTrials;
    z2:.gov.phiInv 1f-1f%nTrials*e;
    (sqrt varSR)*((1f-g)*z1)+g*z2
 };

/ Deflated Sharpe Ratio = PSR with srStar = expected max Sharpe under the null from
/ N trials. DSR >= threshold (.cfg.gov.dsrThreshold, default 0.95) -> the edge
/ survives the multiple-testing penalty. The bar RISES with N (more trials) and
/ FALLS with short n - exactly what kills a tempting small-sample slice.
.gov.deflatedSharpe:{[srHat;n;skew;kurt;nTrials;varSR]
    srStar:.gov.expectedMaxSR[nTrials;varSR];
    .gov.psr[srHat;srStar;n;skew;kurt]
 };

/ Sample moments of a return vector (population sigma, matching .strategy.commodityBT.__perf's
/ `dev`). Returns n, mean, sigma, skew, kurt (non-excess), and the PER-PERIOD Sharpe
/ srPP = mean/sigma (= annualised __perf Sharpe / sqrt(annDays) - the consistency the DSR needs).
.gov.__moments:{[r]
    n:count r;
    if[0=n; :`n`mean`sigma`skew`kurt`srPP!(0;0n;0n;0n;0n;0n)];
    mu:avg r;
    dd:r-mu;
    variance:avg dd*dd;
    sigma:sqrt variance;
    skew:$[sigma>0f; (avg dd*dd*dd)%sigma*sigma*sigma; 0f];
    kurt:$[sigma>0f; (avg dd*dd*dd*dd)%variance*variance; 3f];
    srPP:$[sigma>0f; mu%sigma; 0f];
    `n`mean`sigma`skew`kurt`srPP!(n;mu;sigma;skew;kurt;srPP)
 };

/ ── hypothesis registry + trials ledger (the honest N) ───────────────────────

/ Typed empty tables (so a fresh process / the test suite works with no HDB).
.gov.__emptyHypotheses:{[]
    ([] hypoId:`symbol$(); thesis:(); edgeSource:`symbol$(); instruments:();
        claimedRegimes:(); preRegisteredAt:`timestamp$(); status:`symbol$())
 };
.gov.__emptyTrials:{[]
    ([] trialId:`symbol$(); hypoId:`symbol$(); runAt:`timestamp$(); dataZone:`symbol$();
        dateFrom:`date$(); dateTo:`date$(); nObs:`long$(); regimeAxis:`symbol$();
        regimeBucket:`symbol$(); grossSharpe:`float$(); netSharpe:`float$(); skew:`float$();
        kurtosis:`float$(); maxDD:`float$(); hitRate:`float$(); turnover:`float$();
        paramsJson:(); versionHash:`symbol$())
 };

/ Create the in-memory tables if absent. The on-disk splays are `hypotheses` /
/ `trials`; the in-memory holders are .gov.hypoTbl / .gov.trialTbl (so the query
/ function .gov.trials[hypoId] does not clash with a `trials` global).
.gov.init:{[]
    if[not `hypoTbl in key `.gov; .gov.hypoTbl:.gov.__emptyHypotheses[]];
    if[not `trialTbl in key `.gov; .gov.trialTbl:.gov.__emptyTrials[]];
 };

/ Pre-register a hypothesis (idempotent by hypoId). hypo is a dict with at least
/ `hypoId; optional thesis / edgeSource / instruments / claimedRegimes / status.
/ Returns the hypoId.
.gov.register:{[hypo]
    .gov.init[];
    hid:hypo`hypoId;
    row:`hypoId`thesis`edgeSource`instruments`claimedRegimes`preRegisteredAt`status!(
        hid;
        $[`thesis in key hypo; hypo`thesis; ""];
        $[`edgeSource in key hypo; hypo`edgeSource; `];
        $[`instruments in key hypo; (),hypo`instruments; `symbol$()];
        $[`claimedRegimes in key hypo; (),hypo`claimedRegimes; `symbol$()];
        $[`preRegisteredAt in key hypo; hypo`preRegisteredAt; .z.p];
        $[`status in key hypo; hypo`status; `proposed]);
    .gov.hypoTbl:(delete from .gov.hypoTbl where hypoId=hid) upsert row;
    hid
 };

/ Fetch a registered hypothesis (the row as a dict). Errors if absent.
.gov.hypo:{[hid]
    .gov.init[];
    r:select from .gov.hypoTbl where hypoId=hid;
    if[0=count r; '"gov.hypo: no hypothesis ",string hid];
    first r
 };

/ APPEND one trial to the ledger (NEVER overwrites - the multiple-testing N only
/ grows). trial is a dict; missing fields default. A second logTrial ADDS a row.
.gov.logTrial:{[trial]
    .gov.init[];
    tid:$[`trialId in key trial; trial`trialId; `$"t",string 1+count .gov.trialTbl];
    row:`trialId`hypoId`runAt`dataZone`dateFrom`dateTo`nObs`regimeAxis`regimeBucket`grossSharpe`netSharpe`skew`kurtosis`maxDD`hitRate`turnover`paramsJson`versionHash!(
        tid;
        trial`hypoId;
        $[`runAt in key trial; trial`runAt; .z.p];
        $[`dataZone in key trial; trial`dataZone; `full];
        $[`dateFrom in key trial; trial`dateFrom; 0Nd];
        $[`dateTo in key trial; trial`dateTo; 0Nd];
        $[`nObs in key trial; trial`nObs; 0N];
        $[`regimeAxis in key trial; trial`regimeAxis; `];
        $[`regimeBucket in key trial; trial`regimeBucket; `];
        $[`grossSharpe in key trial; trial`grossSharpe; 0Nf];
        $[`netSharpe in key trial; trial`netSharpe; 0Nf];
        $[`skew in key trial; trial`skew; 0Nf];
        $[`kurtosis in key trial; trial`kurtosis; 0Nf];
        $[`maxDD in key trial; trial`maxDD; 0Nf];
        $[`hitRate in key trial; trial`hitRate; 0Nf];
        $[`turnover in key trial; trial`turnover; 0Nf];
        $[`paramsJson in key trial; trial`paramsJson; "{}"];
        $[`versionHash in key trial; trial`versionHash; `$"v",.qfdm.version]);
    .gov.trialTbl:.gov.trialTbl upsert row;
    tid
 };

/ A family's trials / trial count (= N for deflation).
.gov.trials:{[hid] .gov.init[]; select from .gov.trialTbl where hypoId=hid};
.gov.nTrials:{[hid] count .gov.trials hid};

/ ── HDB persistence (touches real gitignored data; demo / wrapper only) ───────

/ Load the persisted registry + ledger with `get` - NOT \l (no cwd change). Loads
/ the sym domain too. Missing tables -> empty (a fresh HDB is valid).
.gov.open:{[hdbPath]
    symPath:hsym `$hdbPath,"/sym";
    if[0<count key symPath; sym::get symPath];
    hp:hsym `$hdbPath,"/hypotheses/";
    tp:hsym `$hdbPath,"/trials/";
    .gov.hypoTbl:$[0<count key hp; get hp; .gov.__emptyHypotheses[]];
    .gov.trialTbl:$[0<count key tp; get tp; .gov.__emptyTrials[]];
    hdbPath
 };

/ Persist the in-memory registry + ledger back to the splay (enumerate syms via
/ .Q.en). The append-only contract holds across runs: open (load prior) -> logTrial
/ (append in memory) -> flush (write back) loses no prior trial.
.gov.flush:{[hdbPath]
    hp:hsym `$hdbPath;
    (hsym `$hdbPath,"/hypotheses/") set .Q.en[hp;.gov.hypoTbl];
    (hsym `$hdbPath,"/trials/") set .Q.en[hp;.gov.trialTbl];
    hdbPath
 };

/ ── the gate cascade ─────────────────────────────────────────────────────────

/ Build a verdict record (the cascade's output).
.gov.__verdict:{[passed;failedGate;reason;netSharpe;dsr;postHoc;verdict]
    `passedGates`failedGate`reason`netSharpe`dsr`postHoc`verdict!(
        passed;failedGate;reason;netSharpe;dsr;postHoc;verdict)
 };

/ Walk-forward OOS stability: reuse __splits for causal index windows, score each
/ TEST window with __perf (annualised Sharpe). Returns the per-fold OOS Sharpes and
/ the fraction at/above the floor. There is nothing to FIT here (the returns are
/ realised), so the train slice only positions the causal test window - the gate is
/ "is the OOS edge stable across folds", the same causal pattern .alloc.compare uses.
.gov.__walkForward:{[rets;cfg]
    sc:cfg`walkForward;
    nObs:count rets;
    splits:.strategy.commodityBT.__splits[sc`scheme;nObs;sc`trainSpan;sc`testSpan;sc`maxSplits];
    annDays:cfg`annualizationDays;
    floorS:sc`oosSharpeFloor;
    foldSharpes:{[rets;annDays;sp]
        testIdx:(1+sp`trainEndIdx)+til (sp`testEndIdx)-sp`trainEndIdx;
        (.strategy.commodityBT.__perf[rets testIdx;annDays;1f])`sharpe}[rets;annDays] each splits;
    foldSharpes:foldSharpes where not null foldSharpes;
    passFraction:$[0<count foldSharpes; (sum foldSharpes>=floorS)%count foldSharpes; 0f];
    `foldSharpes`passFraction`nFolds!(foldSharpes;passFraction;count foldSharpes)
 };

/ Evaluate the ordered gate cascade for one (hypothesis, net-return vector, family
/ trial count + Sharpe variance). STOPS at the first failure. ev is a dict:
/   hypo (the registered hypothesis dict), rets (NET daily returns of the bucket),
/   bucket (the regime bucket sym), axis (the regime axis sym),
/   nTrials (N from the ledger), varSR (variance of the family's per-period Sharpes),
/   cfg (.cfg.gov or an override).
/ Returns a verdict record: passedGates, failedGate, reason, netSharpe (annualised),
/ dsr, postHoc (bucket not in claimedRegimes -> data-snooped warning), verdict in
/ `reject`/`research`/`regimeConditional`/`pass`.
.gov.evaluate:{[ev]
    hypo:ev`hypo; rets:ev`rets; bucket:ev`bucket;
    nTrials:ev`nTrials; varSR:ev`varSR; cfg:ev`cfg;
    annDays:cfg`annualizationDays;
    claimed:$[`claimedRegimes in key hypo; (),hypo`claimedRegimes; `symbol$()];
    postHoc:not bucket in claimed;
    passed:`symbol$();
    / --- Gate 0: thesis (registered rationale + edgeSource + claimedRegimes) ---
    hasThesis:$[`thesis in key hypo; 0<count hypo`thesis; 0b];
    hasEdge:$[`edgeSource in key hypo; (hypo`edgeSource) in cfg`validEdgeSources; 0b];
    hasClaim:0<count claimed;
    if[not all (hasThesis;hasEdge;hasClaim);
        :.gov.__verdict[passed;`thesis;"missing thesis / valid edgeSource / claimedRegimes";0n;0n;postHoc;`reject]];
    passed,:`thesis;
    / --- Gate 1: cost (net-of-execution annualised Sharpe >= hurdle) ---
    perf:.strategy.commodityBT.__perf[rets;annDays;1f];
    netSharpe:perf`sharpe;
    if[netSharpe < cfg`costHurdleSharpe;
        :.gov.__verdict[passed;`cost;"net Sharpe ",(.gov.__fmt netSharpe)," < hurdle ",.gov.__fmt cfg`costHurdleSharpe;netSharpe;0n;postHoc;`reject]];
    passed,:`cost;
    / --- Gate 2: deflated Sharpe (DSR >= threshold, N from the ledger) ---
    mo:.gov.__moments rets;
    dsr:.gov.deflatedSharpe[mo`srPP;mo`n;mo`skew;mo`kurt;nTrials;varSR];
    if[dsr < cfg`dsrThreshold;
        :.gov.__verdict[passed;`deflatedSharpe;
            "DSR ",(.gov.__fmt dsr)," < ",(.gov.__fmt cfg`dsrThreshold)," (n=",(string mo`n),", N=",(string nTrials),$[postHoc;"; POST-HOC slice";""],")";
            netSharpe;dsr;postHoc;$[postHoc;`regimeConditional;`research]]];
    passed,:`deflatedSharpe;
    / --- Gate 3: walk-forward OOS stability across folds ---
    wf:.gov.__walkForward[rets;cfg];
    wfPass:$[(wf`nFolds)>=cfg[`walkForward;`minFolds]; (wf`passFraction)>=cfg[`walkForward;`minPassFraction]; 0b];
    if[not wfPass;
        :.gov.__verdict[passed;`walkForward;
            "OOS unstable: folds=",(string wf`nFolds)," passFraction=",.gov.__fmt wf`passFraction;
            netSharpe;dsr;postHoc;`research]];
    passed,:`walkForward;
    / --- all gates passed ---
    .gov.__verdict[passed;`none;"all gates passed";netSharpe;dsr;postHoc;`pass]
 };

/ Compact float formatter for reason strings (nulls -> "n/a").
.gov.__fmt:{[x] $[null x; "n/a"; .Q.f[4;x]]};

/ ── the non-optional logging wrapper ─────────────────────────────────────────

/ Group net daily PnL by regime bucket. Returns one dict per bucket:
/ bucket, rets (the bucket's net return vector), dateFrom, dateTo.
.gov.__bucketReturns:{[pnlByDate;regimeLabels;axis]
    labTbl:`date xkey ([] date:regimeLabels`date; bucket:regimeLabels axis);
    m:select from (pnlByDate lj labTbl) where not null bucket;
    buckets:asc distinct m`bucket;
    {[m;b] sub:select from m where bucket=b;
        `bucket`rets`dateFrom`dateTo!(b;sub`pnl;min sub`date;max sub`date)}[m] each buckets
 };

/ THE non-optional path. Logs a trial per regime bucket FIRST (so the ledger's N
/ reflects the multiple-testing count - one trial per slice tried), THEN evaluates
/ each bucket against the now-updated ledger. The backtest ENGINE is untouched -
/ logging is enforced by routing all research through .gov.run, not by an engine hook.
/ Returns one verdict record per bucket (a table). pnlByDate has columns date, pnl
/ (NET daily PnL); regimeLabels is a .regime.series table; axis is e.g. `curveState.
.gov.run:{[hypoId;pnlByDate;regimeLabels;axis]
    .gov.init[];
    cfg:.cfg.gov;
    annDays:cfg`annualizationDays;
    bdata:.gov.__bucketReturns[pnlByDate;regimeLabels;axis];
    / 1. LOG a trial per bucket (append-only -> N grows).
    {[hypoId;axis;annDays;bd]
        r:bd`rets;
        perf:.strategy.commodityBT.__perf[r;annDays;1f];
        mo:.gov.__moments r;
        .gov.logTrial `hypoId`runAt`dataZone`dateFrom`dateTo`nObs`regimeAxis`regimeBucket`grossSharpe`netSharpe`skew`kurtosis`maxDD`hitRate`turnover`paramsJson`versionHash!(
            hypoId;.z.p;`full;bd`dateFrom;bd`dateTo;count r;axis;bd`bucket;
            perf`sharpe;perf`sharpe;mo`skew;mo`kurt;perf`maxDrawdown;perf`hitRate;0n;"{}";`$"v",.qfdm.version)
        }[hypoId;axis;annDays] each bdata;
    / 2. EVALUATE each bucket against the updated ledger (N + V from the family).
    hypo:.gov.hypo hypoId;
    fam:.gov.trials hypoId;
    nTrials:count fam;
    perPeriod:(fam`netSharpe)%sqrt annDays;
    varSR:$[1<count perPeriod; var perPeriod; 0f];
    {[hypo;axis;cfg;nTrials;varSR;bd]
        v:.gov.evaluate `hypo`rets`bucket`axis`nTrials`varSR`cfg!(
            hypo;bd`rets;bd`bucket;axis;nTrials;varSR;cfg);
        `bucket`nObs`verdict`failedGate`netSharpe`dsr`postHoc`reason!(
            bd`bucket;count bd`rets;v`verdict;v`failedGate;v`netSharpe;v`dsr;v`postHoc;v`reason)
        }[hypo;axis;cfg;nTrials;varSR] each bdata
 };

.gov.init[];
-1 "gov.q loaded - .gov.* research governance ready (registry/ledger empty; not opened)";
