/ factor/factor.q - factor decomposition (PCA) of the futures curve (.factor.*) (v0.73, Research OS R8)
/ ----------------------------------------------------------------------------
/ The FIRST new research CAPABILITY built ON the completed spine (R1-R7): PCA of the curve
/ into its principal components (level / slope / curvature) + the residuals. It reads the HDB
/ curve panel DOWNWARD (the same curve accessor the regime layer uses, so the factor view is
/ consistent with the regime view), composes existing functions, and changes no compute path.
/ -
/ q has NO built-in eigendecomposition, so the top-k PCs are computed by DETERMINISTIC power
/ iteration + Hotelling deflation on the covariance matrix. Determinism is PINNED so loadings
/ and scores are byte-identical across runs and known-answer testable: a FIXED initial vector
/ (normalised ones), a FIXED tolerance / max-iter (.cfg.factor), and a SIGN CONVENTION (force a
/ positive loading on the front maturity - eigenvectors are sign-ambiguous, so pinning the sign
/ is mandatory). Registered as a `factor` capability via R2; loads before templates/ (which
/ compose it); never opens the HDB at import.
/ ----------------------------------------------------------------------------

.factor.defaultConfig:{[] .cfg.factor};

/ ── linear-algebra helpers (no built-in eig; power iteration) ────────────────

.factor.__norm:{[v] v % sqrt sum v*v};
.factor.__diag:{[m] m ./: (til count m),'til count m};

/ Top eigenvector/eigenvalue of a (symmetric PSD) matrix by power iteration from a FIXED
/ initial vector, to a FIXED tolerance / max-iter; sign-fixed so the front loading is positive.
.factor.__topEig:{[mat;cfg]
    n:count mat;
    v:.factor.__norm n#1f;
    i:0; done:0b;
    while[(i<cfg`maxIter) and not done;
        w:.factor.__norm mat mmu v;
        done:cfg[`tol] > sqrt sum (w-v)*w-v;
        v:w; i+:1];
    lam:v wsum mat mmu v;
    v:$[v[0]<0f; neg v; v];
    `vec`val!(v;lam)
 };

/ PCA of a T x M panel X (rows = dates, columns = maturities) into the top k PCs.
/ Returns: loadings (M x k - each column a PC), scores (T x k), residuals (T x M = the panel
/ minus its k-factor reconstruction), explainedVar (k-vector, fraction of total variance),
/ colMeans (the demean vector). PURE + deterministic.
.factor.pca:{[X;k;cfg]
    nObs:count X;
    mu:avg X;
    Xc:X -\: mu;
    cmat:(1f%nObs-1) * (flip Xc) mmu Xc;
    trace0:sum .factor.__diag cmat;
    loadings:(); vals:();
    deflated:cmat;
    i:0;
    while[i<k;
        e:.factor.__topEig[deflated;cfg];
        loadings,:enlist e`vec;
        vals,:e`val;
        deflated:deflated - (e`val) * (e`vec) *\: e`vec;
        i+:1];
    loadMat:flip loadings;            / M x k
    scores:Xc mmu loadMat;            / T x k
    residuals:Xc - scores mmu flip loadMat;
    `loadings`scores`residuals`explainedVar`colMeans!(loadMat;scores;residuals;vals%trace0;mu)
 };

/ ── curve panel (reads the HDB curve accessor downward) ──────────────────────

/ A constant-maturity-by-RANK level panel from a curve history (asofDate/tenor/price/...):
/ per date, sort by tenor and take the first M contracts' prices (the same front->deferred rank
/ ordering the regime layer uses). Drops dates with fewer than M live contracts. PURE.
.factor.curvePanel:{[ch;nMat]
    g:0!select prices:price by asofDate from `asofDate`tenor xasc ch;
    g:select from g where nMat<=count each prices;
    `dates`levels!(g`asofDate; nMat sublist/: g`prices)
 };

/ The curve-CHANGE panel: per-maturity price DIFFERENCES across time (drops the seed row).
/ Plain differences (not log returns) - robust to the non-positive prints crude can have (the
/ 2020 negative-WTI settle), and the standard input for curve factor analysis. PURE.
.factor.changePanel:{[levels] 1_ deltas levels};

/ Decompose a commodity's real curve over `dates` (reads the HDB). Convenience wrapper around
/ curvePanel + changePanel + pca; returns the pca dict plus the change-row dates + the levels.
.factor.decompose:{[commodity;dates;cfg]
    ch:.data.hdb.curveHistory[commodity;dates];
    cp:.factor.curvePanel[ch;cfg`nMaturities];
    d:.factor.pca[.factor.changePanel cp`levels; cfg`k; cfg];
    d,`changeDates`levels!(1_ cp`dates; cp`levels)
 };

/ ── register as an R2 `factor` capability ────────────────────────────────────

.factor.contract:`version`requiredIn`requiredOut!(
    1;
    (enlist `panel)!enlist `matrix;
    `loadings`scores`residuals`explainedVar!`matrix`matrix`matrix`floatVec);
.registry.new[`factor; .factor.contract];
.factor.register:{[name;fn;manifest] .registry.register[`factor;name;fn;manifest]};
.factor.get:{[name] .registry.get[`factor;name]};
.factor.list:{[] .registry.list `factor};
.factor.conforms:{[name] .registry.conforms[`factor;name]};

.factor.register[`curvePCA; .factor.pca;
    `contractVersion`description`in`out!(
        1;
        "deterministic power-iteration PCA of the curve-change panel (level/slope/curvature + residuals)";
        (enlist `panel)!enlist `matrix;
        `loadings`scores`residuals`explainedVar!`matrix`matrix`matrix`floatVec)];

-1 "factor.q loaded - .factor.* curve PCA capability ready (registered as `factor/curvePCA)";
