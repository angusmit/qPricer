/ kalmanSchwartzSmith.q - Schwartz-Smith two-factor state-space + Kalman MLE (v0.52)
/ ----------------------------------------------------------------------------
/ Estimates the two-factor commodity model as a linear-Gaussian state space and
/ recovers ALL parameters - crucially the mean-reversion speed kappa AND the
/ volatilities - by Kalman-filter maximum likelihood over a FUTURES PANEL
/ (many dates x tenors). This identifies kappa from the curve DYNAMICS, which a
/ single snapshot (v0.51 calibrateCurve) could not.
/ ----------------------------------------------------------------------------
/ Model (Schwartz-Smith; ln S = chi + xi):
/   d chi = -kappa*chi dt + sigChi dWchi        (short-term deviation, reverts to 0)
/   d xi  = muXi  dt      + sigXi  dWxi          (equilibrium level, BM w/ drift)
/   corr(dWchi,dWxi)=rho
/   ln F(t,T) = e^{-kappa*(T-t)} chi_t + xi_t + A(T-t)
/   A(tau) = muXiStar*tau - (1-e^{-kappa tau})*lamChi/kappa
/            + 0.5*[ (1-e^{-2 kappa tau})*sigChi^2/(2 kappa) + sigXi^2*tau
/                    + 2*(1-e^{-kappa tau})*rho*sigChi*sigXi/kappa ]
/   muXiStar = muXi - lamXi (risk-neutral equilibrium drift). lamChi, lamXi are
/   risk premia; default 0. With lam=0, A(tau) EQUALS schwartz2's deterministic
/   log-futures term (a cross-check, not a dependency) - see the filter test.
/ ----------------------------------------------------------------------------
/ State-space (state alpha=[chi;xi], step dt between observation dates):
/   Transition: alpha_t = c + Tm alpha_{t-1} + eta,  Var(eta)=W
/     Tm=diag(e^{-kappa dt},1); c=(0; muXi*dt)
/     W[0;0]=sigChi^2 (1-e^{-2 kappa dt})/(2 kappa); W[1;1]=sigXi^2 dt
/     W[0;1]=W[1;0]=sigChi sigXi rho (1-e^{-kappa dt})/kappa
/   Measurement: y_t = d_t + Z_t alpha_t + eps,  Var(eps)=V=diag(measSigma^2)
/     y_t = ln(observed futures) for the tenors quoted that date
/     Z_t row for tau: (e^{-kappa tau}, 1);  d_t entry = A(tau)
/   Varying observation dimension per date (missing tenors) is handled naturally
/   by building Z_t/d_t/y_t from only the tenors present that date.
/ ----------------------------------------------------------------------------
/ ASSUMPTIONS (documented): risk premia lamChi=lamXi=0 by default (the synthetic
/ round-trip and example estimate the 6 core params {kappa,muXi,sigChi,sigXi,rho,
/ measSigma}); kappa and the vols are kept in-domain by BOUNDED grid search (no
/ unconstrained transforms needed); the first observation date uses dt=0 (the
/ transition is the identity with zero process noise), so it is a pure update of
/ a diffuse prior. The MLE reuses the repo's grid paradigm (1-D refining grid per
/ coordinate via coordinate descent, scored by the negative log-likelihood and
/ resolved by .calibration.bestCalibrationResult); no new numerical optimiser and
/ no modification to schwartz2/calibrateCurve/parser/optimiser.

.commodity.kalman.__pi:acos neg 1f;
.commodity.kalman.__log2pi:log 2f*.commodity.kalman.__pi;
.commodity.kalman.__priorVar:10f;

/ Identity matrix of dimension n (cycling-take trick).
.commodity.kalman.__ident:{[n] (2#n)#1f,n#0f};

/ Diagonal of a square matrix.
.commodity.kalman.__diag:{[squareMatrix] {[m;i] m[i;i]}[squareMatrix] each til count squareMatrix};

.commodity.kalman.__linspace:{[lo;hi;n]
    $[n<=1; enlist 0.5*lo+hi; lo+(hi-lo)*(til n)%n-1]
 };

.commodity.kalman.defaultParams:{[]
    `kappa`muXi`sigChi`sigXi`correlation`lamChi`lamXi`measSigma!(1.0;0.0;0.30;0.15;0.0;0.0;0.0;0.02)
 };

/ A(tau) deterministic log-futures term, vectorised over a tenor vector.
.commodity.kalman.__Atau:{[taus;params]
    kappa:params`kappa; sigChi:params`sigChi; sigXi:params`sigXi; rho:params`correlation;
    lamChi:params`lamChi; muXiStar:(params`muXi)-params`lamXi;
    eK:exp neg kappa*taus;
    e2K:exp neg 2f*kappa*taus;
    (muXiStar*taus)
        - ((1f-eK)*lamChi%kappa)
        + 0.5*( ((1f-e2K)*sigChi*sigChi%2f*kappa) + (sigXi*sigXi*taus) + (2f*(1f-eK)*rho*sigChi*sigXi%kappa) )
 };

/ Build the state-space matrices for one date (dt to previous date, that date's taus).
.commodity.kalman.buildMatrices:{[params;dt;taus]
    kappa:params`kappa; muXi:params`muXi; sigChi:params`sigChi; sigXi:params`sigXi;
    rho:params`correlation; measSig:params`measSigma;
    eKdt:exp neg kappa*dt;
    Tm:(eKdt,0f;0f,1f);
    cVec:(0f;muXi*dt);
    W00:$[dt>0f; (sigChi*sigChi)*(1f-exp neg 2f*kappa*dt)%2f*kappa; 0f];
    W11:(sigXi*sigXi)*dt;
    W01:$[dt>0f; (sigChi*sigXi*rho)*(1f-eKdt)%kappa; 0f];
    Wm:(W00,W01;W01,W11);
    nObs:count taus;
    Zm:flip (exp neg kappa*taus;nObs#1f);
    dVec:.commodity.kalman.__Atau[taus;params];
    Vm:(measSig*measSig)*.commodity.kalman.__ident nObs;
    `Tm`c`W`Z`d`V!(Tm;cVec;Wm;Zm;dVec;Vm)
 };

/ One Kalman predict+update step on a pre-built matrix bundle (generic; works for
/ any state/obs dimension). Accumulates the Gaussian log-likelihood.
.commodity.kalman.__kalmanStep:{[state;bundle]
    a:state`a; P:state`P; ll:state`loglik;
    Tm:bundle`Tm; cVec:bundle`c; Wm:bundle`W; Zm:bundle`Z; dVec:bundle`d; Vm:bundle`V; y:bundle`y;
    aPred:cVec+Tm mmu a;
    PPred:(Tm mmu P mmu flip Tm)+Wm;
    yhat:dVec+Zm mmu aPred;
    vInnov:y-yhat;
    Sm:(Zm mmu PPred mmu flip Zm)+Vm;
    SmInv:inv Sm;
    cholSm:.correlation.__cholesky Sm;
    logDet:2f*sum log .commodity.kalman.__diag cholSm;
    nObs:count y;
    quad:sum vInnov*SmInv mmu vInnov;
    llStep:neg 0.5*((nObs*.commodity.kalman.__log2pi)+logDet+quad);
    Kgain:PPred mmu (flip Zm) mmu SmInv;
    aNew:aPred+Kgain mmu vInnov;
    PNew:PPred-Kgain mmu Zm mmu PPred;
    `a`P`loglik!(aNew;PNew;ll+llStep)
 };

/ Run the filter over a list of matrix bundles from a seed (a0,P0). Returns total
/ log-likelihood and the per-step filtered states.
.commodity.kalman.filterBundles:{[bundles;a0;P0]
    seed:`a`P`loglik!(a0;P0;0f);
    states:.commodity.kalman.__kalmanStep\[seed;bundles];
    `loglik`states!((last states)`loglik;states)
 };

.commodity.kalman.__validatePanel:{[panel]
    if[not 98h=type panel; '"kalman: panel must be a table"];
    if[not all `obsDate`tau`logF in cols panel; '"kalman: panel needs obsDate, tau, logF columns"];
    if[2>count distinct panel`obsDate; '"kalman: need >=2 observation dates"];
 };

/ Filter the Schwartz-Smith model over a (obsDate, tau, logF) panel.
.commodity.kalman.filter:{[panel;params]
    .commodity.kalman.__validatePanel panel;
    dates:asc distinct panel`obsDate;
    nDates:count dates;
    dtVec:0f,(`float$1_deltas dates)%365f;
    bundles:();
    di:0;
    while[di<nDates;
        d:dates di;
        sub:`tau xasc select tau,logF from panel where obsDate=d;
        mats:.commodity.kalman.buildMatrices[params;dtVec di;sub`tau];
        bundles,:enlist mats,enlist[`y]!enlist sub`logF;
        di+:1];
    firstLogF:exec logF from panel where obsDate=first dates;
    a0:(0f;avg firstLogF);
    pv:.commodity.kalman.__priorVar;
    P0:(pv,0f;0f,pv);
    fb:.commodity.kalman.filterBundles[bundles;a0;P0];
    states:fb`states;
    `loglik`dates`chi`xi`status!(fb`loglik;dates;states[;`a][;0];states[;`a][;1];`OK)
 };

/ Build a Kalman panel (obsDate, tau, logF) from a .parser.crude.curveHistory
/ table, excluding non-positive prices (log domain).
.commodity.kalman.panelFromCurveHistory:{[curveHistory]
    if[not 98h=type curveHistory; '"kalman panel: curveHistory must be a table"];
    if[not all `asofDate`tenor`price in cols curveHistory; '"kalman panel: needs asofDate, tenor, price"];
    positiveRows:select from curveHistory where price>0f;
    if[0=count positiveRows; '"kalman panel: no positive prices"];
    `obsDate xasc select obsDate:asofDate, tau:`float$tenor, logF:log `float$price from positiveRows
 };

/ ----------------------------------------------------------------------------
/ Maximum-likelihood estimation (coordinate descent of negative log-likelihood)
/ ----------------------------------------------------------------------------

.commodity.kalman.defaultEstCfg:{[]
    `bounds`init`gridSteps`refineRounds`refineShrink`nSweeps!(
        `kappa`muXi`sigChi`sigXi`correlation`measSigma!((0.05 5.0);(-0.5 0.5);(0.01 1.5);(0.01 1.0);(-0.95 0.95);(0.0001 0.3));
        `kappa`muXi`sigChi`sigXi`correlation`measSigma!(1.0;0.0;0.30;0.15;0.0;0.02);
        11;
        4;
        0.4;
        5)
 };

.commodity.kalman.__freeNames:`kappa`muXi`sigChi`sigXi`correlation`measSigma;

/ Build a full params dict from a free-parameter vector (lamChi=lamXi=0 fixed).
.commodity.kalman.__vecToParams:{[vec]
    (.commodity.kalman.defaultParams[]),.commodity.kalman.__freeNames!vec
 };

/ Negative log-likelihood of a free-parameter vector (errors -> +inf so rejected).
.commodity.kalman.__negLogLik:{[panel;vec]
    params:.commodity.kalman.__vecToParams vec;
    ll:@[{[pnl;prm] (.commodity.kalman.filter[pnl;prm])`loglik}[panel;];params;{[e] -0w}];
    $[(null ll)|ll=0w; 0w; neg ll]
 };

/ 1-D refining grid minimiser of g over [lo,hi]; reuses bestCalibrationResult.
.commodity.kalman.__refine1D:{[g;lo;hi;gridSteps;refineRounds;shrink]
    loHi:(lo;hi);
    bestX:0n; bestVal:0w; totalIter:0;
    roundIdx:0;
    while[roundIdx<refineRounds;
        grid:.commodity.kalman.__linspace . loHi,gridSteps;
        vals:{[gg;x] .[gg;enlist x;{[e] 0w}]}[g;]'[grid];
        candTable:([] cand:grid; rmse:vals; status:(count grid)#`OK);
        bestRow:.calibration.bestCalibrationResult candTable;
        totalIter+:count grid;
        if[(bestRow`rmse)<bestVal; bestVal:bestRow`rmse; bestX:bestRow`cand];
        half:0.5*shrink*loHi[1]-loHi 0;
        loHi:((lo|bestX-half);hi&bestX+half);
        roundIdx+:1];
    `bestX`bestVal`iterations!(bestX;bestVal;totalIter)
 };

.commodity.kalman.estimate:{[panel;calCfg]
    .commodity.kalman.__validatePanel panel;
    cfg:.commodity.kalman.defaultEstCfg[];
    if[count calCfg; cfg:cfg,calCfg];
    freeNames:.commodity.kalman.__freeNames;
    bounds:cfg`bounds;
    current:`float$cfg[`init] freeNames;
    nSweeps:cfg`nSweeps;
    negLL:.commodity.kalman.__negLogLik[panel;];
    totalIter:0;
    sweepIdx:0;
    while[sweepIdx<nSweeps;
        ci:0;
        while[ci<count freeNames;
            loHi:bounds freeNames ci;
            g:{[fn;cur;idx;x] fn @[cur;idx;:;x]}[negLL;current;ci;];
            best:.commodity.kalman.__refine1D[g;loHi 0;loHi 1;cfg`gridSteps;cfg`refineRounds;cfg`refineShrink];
            current[ci]:best`bestX;
            totalIter+:best`iterations;
            ci+:1];
        sweepIdx+:1];
    estParams:.commodity.kalman.__vecToParams current;
    finalFilter:.commodity.kalman.filter[panel;estParams];
    `estimatedParams`loglik`chi`xi`dates`iterations`status!(
        estParams;finalFilter`loglik;finalFilter`chi;finalFilter`xi;finalFilter`dates;totalIter;`OK)
 };

/ ----------------------------------------------------------------------------
/ Synthetic panel simulator (for round-trip testing / demonstration)
/ ----------------------------------------------------------------------------
/ Simulate chi/xi forward over nDates dates spaced dayStep calendar days apart
/ (dt=dayStep/365, matching what filter[] derives), observe ln-futures at the
/ given (fixed) tenors with measurement noise. Returns the (obsDate,tau,logF)
/ panel plus the true chi/xi paths.
.commodity.kalman.simulatePanel:{[params;dayStep;nDates;taus;startDate;xi0;seed]
    if[nDates<2; '"simulatePanel: need >=2 dates"];
    if[dayStep<=0; '"simulatePanel: dayStep must be positive"];
    dt:dayStep%365f;
    mats:.commodity.kalman.buildMatrices[params;dt;taus];
    cholW:.correlation.__cholesky mats`W;
    kappa:params`kappa; muXi:params`muXi; measSig:params`measSigma;
    eKdt:exp neg kappa*dt;
    Atau:.commodity.kalman.__Atau[taus;params];
    loadings:exp neg kappa*taus;
    nObs:count taus;
    etaNormals:.commodity.kalman.__matrixReshape[.montecarlo.__generateNormals[2*nDates;seed];nDates;2];
    epsNormals:.commodity.kalman.__matrixReshape[.montecarlo.__generateNormals[nObs*nDates;seed+104729];nDates;nObs];
    dates:startDate+dayStep*til nDates;
    chi:0f; xi:xi0;
    chiPath:(); xiPath:(); panelRows:();
    ti:0;
    while[ti<nDates;
        if[ti>0;
            etaPair:cholW mmu etaNormals ti;
            chi:(eKdt*chi)+etaPair 0;
            xi:(xi+muXi*dt)+etaPair 1];
        logFvec:(loadings*chi)+xi+Atau+measSig*epsNormals ti;
        panelRows,:enlist `obsDate`tau`logF!((nObs#dates ti);taus;logFvec);
        chiPath,:chi; xiPath,:xi;
        ti+:1];
    panel:`obsDate xasc raze {[r] flip r}each panelRows;
    `panel`chi`xi`dates!(panel;chiPath;xiPath;dates)
 };

/ Reshape a flat vector into an r x c matrix (row-major).
.commodity.kalman.__matrixReshape:{[flatVec;r;c] (r;c)#flatVec};
