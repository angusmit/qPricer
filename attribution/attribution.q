/ attribution/attribution.q - PnL explain + bucketed curve risk (.attribution.*) (v0.80, Research OS R15)
/ ----------------------------------------------------------------------------
/ ARCHITECTURE.md Part II 11.8(g) / 13-R15. Turns "explain why this makes money" from prose into a
/ QUANTITATIVE attribution: given a replay run record (R12) + the curve (R10), decompose the realized
/ PnL into commodity-market components - level / slope / curvature / carry / residual - and compute the
/ bucketed curve risk on the positions. This is the quantitative form of the "name which of the three
/ edges" discipline (Sec 10): mostly CARRY = a harvested risk premium; mostly SLOPE/curvature = a
/ structural relative-value normalisation; a large RESIDUAL = unexplained (a red flag for noise/overfit).
/ -
/ LAYERING (the correction): the attribution is a NEW HIGH layer ABOVE backtest/ - it CONSUMES a replay
/ run record (which lives in backtest/), so it canNOT live in analytics/ (analytics/ is BELOW backtest/;
/ that would be an illegal upward dependency). It reads backtest/ (the run record - a data structure),
/ curve/ (R10's curve + the parallel/slope/butterfly shock operators - the decomposition BASIS - and
/ rollYield), and state/ (R9) - all downward. Loads after the backtest block; opens nothing at import.
/ ADDITIVE: it does NOT touch the existing analytics/ risk functions and does NOT rewire cards/.
/ -
/ THE DECOMPOSITION (reconciles by construction). Per step the active contract's price change splits as
/ P_t(tau_t) - P_{t-1}(tau_{t-1}) = [P_t(tau_t) - P_t(tau_{t-1})]  (the CARRY / roll-down along the
/ current curve as the contract ages) + [P_t(tau_{t-1}) - P_{t-1}(tau_{t-1})] (the curve SHIFT at the
/ fixed tenor). The shift is projected onto R10's parallel/slope/butterfly shapes (normal equations) ->
/ level/slope/curvature; whatever the basis does not span is the curve RESIDUAL. residual is DEFINED as
/ realizedStepPnL - (level+slope+curvature+carry), so it also absorbs costs/financing - and the five
/ ALWAYS reconcile to the realized total (like R13's pnlTies). A large residual FRACTION is the
/ quantitative "name which edge, or it's noise."
/ -
/ RESERVED-NAME NOTE: `asOf` NEVER `asof`; `comm` NOT `commodity`; builtins not shadowed (lsq/mmu/inv/
/ sum/wsum/avg/dev/var/first/last/bin/binr). The basis is derived FROM R10's shock operators (not reinvented).
/ ----------------------------------------------------------------------------

.attribution.defaultConfig:{[] .cfg.attribution};

/ ── helpers ──────────────────────────────────────────────────────────────────

/ Linear interpolation of a curve (sorted ascending tenors xs -> prices ys) at tenor x (flat outside).
.attribution.__interp:{[xs;ys;x]
    $[x<=first xs; first ys;
      x>=last xs; last ys;
      [j:xs bin x; x0:xs j; x1:xs j+1; (ys j)+((ys j+1)-ys j)*(x-x0)%x1-x0]]
 };

/ The level/slope/curvature shape vectors on a curve, derived FROM R10's shock operators (a unit shock
/ minus the base price): parallel -> 1s, slope -> tenor-first, butterfly -> the 4x(1-x) bow.
.attribution.__shapes:{[curve]
    px:`float$curve`price;
    (((.curve.shock.parallel[curve;1f])`price)-px;
     ((.curve.shock.slope[curve;1f])`price)-px;
     ((.curve.shock.butterfly[curve;1f])`price)-px)
 };

/ Least-squares coefficients (3-vec) for dP ~ B*coefs (B is M x 3) via the normal equations; 0 0 0f on
/ too-few rows / a singular system (then the whole shift lands in the residual).
.attribution.__solve:{[B;dP]
    $[3>count dP; 0 0 0f;
        @[{[B;dP] bt:flip B; (inv bt mmu B) mmu bt mmu dP}[B;]; dP; {[e] 0 0 0f}]]
 };

/ Per-step component tuple (level;slope;curvature;carry) in the run's return-PnL units. cc is the
/ date->curve-object cache; t is the step index (t=0 -> all zero: no prior curve).
.attribution.__stepComp:{[cc;steps;prevPos;t]
    if[t=0; :0 0 0 0f];
    d:steps[t]`stepDate; p:steps[t-1]`stepDate; c:steps[t]`activeContract;
    ccur:cc d; cprev:cc p;
    if[((::)~ccur) or (::)~cprev; :0 0 0 0f];
    curC:ccur`curve; curP:cprev`curve;
    cymP:curP`contractYM; cymC:curC`contractYM;
    iP:cymP?c; iC:cymC?c;
    if[(iP=count cymP) or iC=count cymC; :0 0 0 0f];     / active not on both curves -> all residual
    Pp:`float$(curP`price) iP;
    if[Pp<=0f; :0 0 0 0f];
    Pt:`float$(curC`price) iC;
    tauP:`float$(curP`tenor) iP;
    shp:.attribution.__shapes curP; L:shp 0; S:shp 1; Cv:shp 2;
    / curve SHIFT at FIXED tenor: interpolate the CURRENT curve at each prior tenor, minus the prior
    / price (so the roll-down/aging is NOT in the projection - that is the carry term, computed below).
    ccurT:`float$curC`tenor; ccurP:`float$curC`price;
    tenP:`float$curP`tenor;
    dP:(.attribution.__interp[ccurT;ccurP;] each tenP)-`float$curP`price;
    coefs:.attribution.__solve[flip (L;S;Cv); dP];
    / carry / roll-down: the active contract's actual new price vs the current curve at its OLD tenor.
    carryMove:Pt - .attribution.__interp[ccurT; ccurP; tauP];
    pp:prevPos t;
    ((pp*(coefs 0)*(L iP)%Pp);                            / level
     (pp*(coefs 1)*(S iP)%Pp);                            / slope
     (pp*(coefs 2)*(Cv iP)%Pp);                           / curvature
     (pp*carryMove%Pp))                                   / carry (roll-down)
 };

/ ── PnL decomposition (.attribution.pnl) ─────────────────────────────────────

/ Decompose a replay run's realized PnL into level/slope/curvature/carry/residual (+ the residual
/ fraction). The five reconcile to the run's realized total by construction (residual is the plug).
.attribution.pnl:{[run]
    steps:run`steps; comm:(run`meta)`commodity;
    n:count steps;
    if[0=n; '"attribution.pnl: empty run record"];
    udts:asc distinct steps`stepDate;
    cc:udts!{[comm;dd] @[{.curve.build[x;y]}[;comm];dd;{[e] (::)}]}[comm] each udts;
    prevPos:0f^prev steps`position;                       / the position HELD during each step (pre-fill)
    comps:.attribution.__stepComp[cc;steps;prevPos;] each til n;
    lev:comps[;0]; slo:comps[;1]; cur:comps[;2]; car:comps[;3];
    sp:steps`stepPnl;
    resid:sp-(lev+slo+cur+car);
    total:sum sp;
    `level`slope`curvature`carry`residual`total`residualFraction`reconciled!(
        sum lev; sum slo; sum cur; sum car; sum resid; total;
        $[0f=abs total; 0n; (abs sum resid)%abs total];
        (.cfg.attribution`reconcileTol)>abs total-(sum lev)+(sum slo)+(sum cur)+(sum car)+sum resid)
 };

/ ── bucketed curve risk (.attribution.risk) ──────────────────────────────────

/ The bucketed curve delta (the end position's $ sensitivity to a unit move per tenor bucket) + the
/ roll-down exposure (position x the end curve's rollYield) + the calendar-spread exposure on the
/ end-of-run position. (Inter-commodity basis exposure is out of scope for a single-commodity run.)
.attribution.risk:{[run]
    steps:run`steps; comm:(run`meta)`commodity;
    buckets:.cfg.attribution`buckets;
    endDate:max steps`stepDate;
    endPos:`float$last steps`position;
    endActive:last steps`activeContract;
    cv:.curve.build[endDate;comm]; curC:cv`curve;
    iA:(curC`contractYM)?endActive;
    tauA:$[iA<count curC`contractYM; `float$(curC`tenor) iA; 0n];
    bucketDelta:(`$"b",/:string til count buckets)!(count buckets)#0f;
    if[not null tauA; bucketDelta[(`$"b",string buckets bin tauA)]:endPos];
    `endDate`endContract`endPosition`tenor`buckets`bucketDelta`rollDownExposure`calendarSpreadExposure!(
        endDate; endActive; endPos; tauA; buckets; bucketDelta;
        endPos*(cv`features)`rollYield;                   / roll-down (carry) exposure
        endPos)                                           / single front leg: net spread exposure = the position
 };

/ ── register as an R2 `attribution capability ────────────────────────────────

.attribution.contract:`version`requiredIn`requiredOut!(
    1;
    `run`commodity!`dict`symbol;
    `level`slope`curvature`carry`residual!`float`float`float`float`float);
.registry.new[`attribution; .attribution.contract];
.attribution.register:{[name;fn;manifest] .registry.register[`attribution;name;fn;manifest]};
.attribution.get:{[name] .registry.get[`attribution;name]};
.attribution.list:{[] .registry.list `attribution};
.attribution.conforms:{[name] .registry.conforms[`attribution;name]};

.attribution.register[`pnlAttribution; .attribution.pnl;
    `contractVersion`description`in`out!(
        1;
        "decompose a replay run's realized PnL into level/slope/curvature/carry/residual (R10 shock basis) + the residual fraction; reconciles to the realized total";
        `run`commodity!`dict`symbol;
        `level`slope`curvature`carry`residual!`float`float`float`float`float)];

-1 "attribution.q loaded - .attribution.* PnL decomposition + bucketed curve risk ready (a HIGH layer above backtest/; not opened)";
