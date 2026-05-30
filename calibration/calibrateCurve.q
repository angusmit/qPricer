/ calibrateCurve.q - fit a commodity term-structure model to a market forward
/ curve in the (tenor, price) shape produced by .parser.crude.curveAt (v0.51)
/ ----------------------------------------------------------------------------
/ Public:  .commodity.calibrateCurve[marketCurve; modelType; calCfg]
/ Helpers: .commodity.curveCal.*  (private __ helpers + defaultCalCfg)
/ ----------------------------------------------------------------------------
/ REUSE, NOT REWRITE: the repository's calibration paradigm (Heston/Bates/SABR)
/ is a GRID SEARCH scored by .objective.* and resolved by
/ .calibration.bestCalibrationResult; there is no gradient optimiser. This
/ module reuses exactly that - a refining grid scored by .objective.rmse with
/ .calibration.bestCalibrationResult picking the best row - and CALLS the
/ existing .commodity.schwartz2.futuresCurve (vectorised over the tenor vector).
/ It never modifies the model, optimiser, or objective modules.
/ ----------------------------------------------------------------------------
/ MODEL & IDENTIFIABILITY (schwartz2, the v0.51 target - for crude):
/   log F(tau) = X0*exp(-kappa*tau) + Y0 + muY*tau + 0.5*v(tau)
/ For a FIXED kappa this is LINEAR in (X0, Y0, muY): the columns exp(-kappa*tau),
/ 1, tau span the deterministic part and 0.5*v(tau) is a known convexity offset
/ (extracted from the existing futures function as log F(X0=0,Y0=0,muY=0)). So we
/ solve (X0, Y0, muY) by EXACT linear least squares (normal equations via q inv)
/ for each candidate kappa, and search only the single nonlinear parameter kappa
/ on a refining 1-D grid. This avoids the shallow X0/kappa/muY trade-off that
/ defeats an axis-aligned 3-D grid and recovers known params cleanly.
/   FREE  : shortFactor0 (X0, initial short factor / convenience-yield deviation),
/           longFactor0 (Y0, long-run log-price level), longDrift (muY), and
/           meanReversionSpeed (kappa, the 1-D searched parameter).
/   FIXED : shortVolatility, longVolatility, correlation (at calCfg values) - a
/           single curve snapshot under-identifies the vols (they enter only via
/           the small convexity term), so they are held at calCfg values.
/ modelType `schwartz / `mrjump are accepted by the signature but reserved
/ (single-curve one-factor / jump fits are out of v0.51 scope) and raise a
/ controlled error.
/ ----------------------------------------------------------------------------
/ NEGATIVE-PRICE GUARD: the Schwartz models are log-price; a curve with any
/ price <= 0 (e.g. the April-2020 WTI negative-settle regime) is out of domain
/ and raises a controlled error rather than being fed to the lognormal fit.
/ ----------------------------------------------------------------------------
/ ECONOMIC READ: from the fitted curve we report impliedFrontSlope =
/ d logF/dtau at the front and netConvenienceYield = riskFreeRate - frontSlope.
/ Backwardation (front > deferred, slope < 0) => convenience yield > carry rate.
/ At least 3 curve points are required (3 free linear parameters per kappa).

.commodity.curveCal.defaultCalCfg:{[]
    .cfg.calib.curve
 };

/ Validate the market curve and return it sorted ascending by tenor (unkeyed dict
/ of tenor/price float vectors).
.commodity.curveCal.__validateCurve:{[marketCurve]
    if[not 98h=type marketCurve; '"calibrateCurve: marketCurve must be a table"];
    if[3>count marketCurve; '"calibrateCurve: need at least 3 curve points (3-parameter linear fit per kappa)"];
    if[not all `tenor`price in cols marketCurve; '"calibrateCurve: marketCurve needs tenor and price columns"];
    sorted:`tenor xasc marketCurve;
    tenors:`float$sorted`tenor;
    prices:`float$sorted`price;
    if[any prices<=0f; '"calibrateCurve: non-positive price in curve - out of lognormal domain (e.g. negative settles)"];
    if[any tenors<=0f; '"calibrateCurve: non-positive tenor in curve"];
    `tenor`price!(tenors;prices)
 };

/ Build a schwartz2 params dict from a kappa, a muY, and the fixed vols.
.commodity.curveCal.__params:{[kappa;muY;fixedCfg]
    `meanReversionSpeed`shortVolatility`longVolatility`longDrift`correlation!(
        kappa;fixedCfg`shortVolatility;fixedCfg`longVolatility;muY;fixedCfg`correlation)
 };

/ Exact linear least-squares fit of (X0, Y0, muY) for a fixed kappa, then the
/ resulting model curve and rmse. ctx carries fixedCfg, tenorVec, marketPrice,
/ logMarket. Returns dict rmse, shortFactor0, longFactor0, longDrift, modelCurve.
.commodity.curveCal.__fitAtKappa:{[kappa;ctx]
    fixedCfg:ctx`fixedCfg;
    tenorVec:ctx`tenorVec;
    marketPrice:ctx`marketPrice;
    / Convexity 0.5*v(tau) = log F(X0=0,Y0=0,muY=0) from the existing futures fn.
    convexityParams:.commodity.curveCal.__params[kappa;0f;fixedCfg];
    convexity:log .commodity.schwartz2.futuresCurve[0f;0f;convexityParams;tenorVec];
    yTarget:(ctx`logMarket)-convexity;
    designMatrix:flip (exp neg kappa*tenorVec;(count tenorVec)#1f;tenorVec);
    designT:flip designMatrix;
    beta:(inv designT mmu designMatrix) mmu designT mmu yTarget;
    x0:beta 0; y0:beta 1; muY:beta 2;
    fitParams:.commodity.curveCal.__params[kappa;muY;fixedCfg];
    modelCurve:.commodity.schwartz2.futuresCurve[x0;y0;fitParams;tenorVec];
    `rmse`shortFactor0`longFactor0`longDrift`modelCurve!(
        .objective.rmse[modelCurve;marketPrice];x0;y0;muY;modelCurve)
 };

.commodity.curveCal.__linspace:{[lo;hi;n]
    $[n<=1; enlist 0.5*lo+hi; lo+(hi-lo)*(til n)%n-1]
 };

/ Refining 1-D grid search over kappa. evalRmse is unary kappa -> rmse. Reuses
/ .calibration.bestCalibrationResult to pick the min-rmse kappa each round.
.commodity.curveCal.__searchKappa:{[evalRmse;kappaRange;gridSteps;refineRounds;shrink;kappaFloor]
    loHi:kappaRange;
    bestKappa:0n;
    bestRmse:0w;
    totalIter:0;
    roundIdx:0;
    while[roundIdx<refineRounds;
        kappaGrid:.commodity.curveCal.__linspace . loHi,gridSteps;
        rmseVec:{[ev;k] .[ev;enlist k;{[e] 0w}]}[evalRmse;]'[kappaGrid];
        candTable:([] meanReversionSpeed:kappaGrid; rmse:rmseVec; status:(count kappaGrid)#`OK);
        bestRow:.calibration.bestCalibrationResult candTable;
        totalIter+:count kappaGrid;
        if[(bestRow`rmse)<bestRmse;
            bestRmse:bestRow`rmse;
            bestKappa:bestRow`meanReversionSpeed];
        halfWidth:0.5*shrink*loHi[1]-loHi 0;
        loHi:((kappaFloor|kappaRange[0]|bestKappa-halfWidth);kappaRange[1]&bestKappa+halfWidth);
        roundIdx+:1];
    `bestKappa`bestRmse`iterations!(bestKappa;bestRmse;totalIter)
 };

.commodity.calibrateCurve:{[marketCurve;modelType;calCfg]
    cleanCurve:.commodity.curveCal.__validateCurve marketCurve;
    if[not modelType in `schwartz2`schwartz`mrjump;
        '"calibrateCurve: unknown modelType ",string modelType];
    if[not modelType=`schwartz2;
        '"calibrateCurve: modelType ",(string modelType)," reserved; v0.51 calibrates schwartz2 only"];
    cfg:.commodity.curveCal.defaultCalCfg[];
    if[count calCfg; cfg:cfg,calCfg];
    tenorVec:cleanCurve`tenor;
    marketPrice:cleanCurve`price;
    fixedCfg:`shortVolatility`longVolatility`correlation#cfg;
    ctx:`fixedCfg`tenorVec`marketPrice`logMarket!(fixedCfg;tenorVec;marketPrice;log marketPrice);
    evalRmse:{[ctxL;kappa] (.commodity.curveCal.__fitAtKappa[kappa;ctxL])`rmse}[ctx;];
    searchResult:.commodity.curveCal.__searchKappa[evalRmse;cfg`meanReversionSpeedRange;cfg`gridSteps;cfg`refineRounds;cfg`refineShrink;cfg`kappaFloor];
    bestKappa:searchResult`bestKappa;
    fit:.commodity.curveCal.__fitAtKappa[bestKappa;ctx];
    fittedPrices:fit`modelCurve;
    fitRmse:fit`rmse;
    calibratedParams:`shortFactor0`longFactor0`meanReversionSpeed`shortVolatility`longVolatility`longDrift`correlation!(
        fit`shortFactor0;fit`longFactor0;bestKappa;fixedCfg`shortVolatility;fixedCfg`longVolatility;fit`longDrift;fixedCfg`correlation);
    perTenorError:([] tenor:tenorVec; marketPrice:marketPrice; modelPrice:fittedPrices; error:fittedPrices-marketPrice);
    fittedCurve:([] tenor:tenorVec; modelPrice:fittedPrices);
    riskFreeRate:cfg`riskFreeRate;
    impliedFrontSlope:((log fittedPrices 1)-log fittedPrices 0)%tenorVec[1]-tenorVec 0;
    netConvenienceYield:riskFreeRate-impliedFrontSlope;
    `calibratedParams`fitRmse`perTenorError`fittedCurve`iterations`modelType`riskFreeRate`impliedFrontSlope`netConvenienceYield`status!(
        calibratedParams;fitRmse;perTenorError;fittedCurve;searchResult`iterations;modelType;riskFreeRate;impliedFrontSlope;netConvenienceYield;`OK)
 };

/ ----------------------------------------------------------------------------
/ CONVENIENCE-YIELD TIME SERIES (v0.52, Part A)
/ ----------------------------------------------------------------------------
/ Calibrate the curve per as-of date across a .parser.crude.curveHistory panel
/ to produce a convenience-yield time series. kappa is FIXED at calCfg.kappa: a
/ single curve snapshot does not identify the mean-reversion speed (freeing it
/ per date injects noise), so we hold it constant and free only the log-linear
/ params (X0, Y0, muY) -- exactly the v0.51 fit with a degenerate kappa range
/ (meanReversionSpeedRange=(kappa;kappa)), reusing .commodity.calibrateCurve
/ unchanged. Identifying kappa from the curve DYNAMICS is Part B's job (Kalman
/ MLE); the example feeds that kappa back here for a cleaner series.
/ Dates whose curve contains any non-positive price are SKIPPED (this excludes
/ the April-2020 WTI negative-settle window automatically); per-date fit errors
/ (e.g. fewer than 3 contracts) are isolated and also recorded as skipped.
/ Returns a dict: series (time-series table), skipped (date list), transitions
/ (regime-flip table backwardation<->contango).

.commodity.curveCal.defaultSeriesCfg:{[]
    .cfg.calib.curveSeries
 };

.commodity.curveCal.__seriesEmpty:{[]
    ([] asofDate:`date$(); frontPrice:`float$(); netConvenienceYield:`float$(); slope:`float$();
        regime:`symbol$(); fitRmse:`float$(); X0:`float$(); Y0:`float$(); muY:`float$(); status:`symbol$())
 };

/ Fit one date's curve; returns a status-tagged dict (OK row or skipped marker).
/ Param is asofVal (NOT asofDate) to avoid colliding with the qSQL column name.
.commodity.curveCal.__fitOneDate:{[asofVal;curveHistory;perDateCfg;rateVal]
    sub:`tenor xasc select tenor,price from curveHistory where asofDate=asofVal;
    prices:sub`price;
    if[any prices<=0f; :`status`asofDate!(`skipped;asofVal)];
    res:@[.commodity.calibrateCurve[sub;`schwartz2;];perDateCfg;{[e] `status`errorMessage!(`ERROR;e)}];
    if[`ERROR~res`status; :`status`asofDate!(`skipped;asofVal)];
    cp:res`calibratedParams;
    cyVal:res`netConvenienceYield;
    `status`asofDate`frontPrice`netConvenienceYield`slope`regime`fitRmse`X0`Y0`muY!(
        `OK;asofVal;first prices;cyVal;res`impliedFrontSlope;
        $[cyVal>rateVal;`backwardation;`contango];res`fitRmse;cp`shortFactor0;cp`longFactor0;cp`longDrift)
 };

.commodity.curveCal.convenienceYieldSeries:{[curveHistory;calCfg]
    if[not 98h=type curveHistory; '"convenienceYieldSeries: curveHistory must be a table"];
    if[not all `asofDate`tenor`price in cols curveHistory; '"convenienceYieldSeries: needs asofDate, tenor, price columns"];
    if[0=count curveHistory; '"convenienceYieldSeries: empty curveHistory"];
    cfg:.commodity.curveCal.defaultSeriesCfg[];
    if[count calCfg; cfg:cfg,calCfg];
    kappaFixed:cfg`kappa;
    rateVal:cfg`riskFreeRate;
    / Per-date calCfg: fix kappa via a degenerate search range, single evaluation.
    perDateCfg:(`shortVolatility`longVolatility`correlation`riskFreeRate#cfg),
        `meanReversionSpeedRange`gridSteps`refineRounds!((kappaFixed;kappaFixed);1;1);
    dates:asc distinct curveHistory`asofDate;
    results:.commodity.curveCal.__fitOneDate[;curveHistory;perDateCfg;rateVal] each dates;
    statuses:results[;`status];
    okResults:results where `OK=statuses;
    skippedDates:results[;`asofDate] where `skipped=statuses;
    seriesCols:`asofDate`frontPrice`netConvenienceYield`slope`regime`fitRmse`X0`Y0`muY`status;
    seriesTable:$[count okResults;
        flip seriesCols!{[colName;rowDicts] rowDicts[;colName]}[;okResults] each seriesCols;
        .commodity.curveCal.__seriesEmpty[]];
    seriesTable:`asofDate xasc seriesTable;
    / Regime transitions: dates where the regime flips vs the previous date.
    transitions:([] asofDate:`date$(); fromRegime:`symbol$(); toRegime:`symbol$());
    if[1<count seriesTable;
        regimes:seriesTable`regime;
        dts:seriesTable`asofDate;
        flipIdx:(where regimes<>prev regimes) except 0;
        transitions:([] asofDate:dts flipIdx; fromRegime:regimes flipIdx-1; toRegime:regimes flipIdx)];
    `series`skipped`transitions!(seriesTable;skippedDates;transitions)
 };
