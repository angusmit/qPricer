/ curve/curve.q - the curve engine (.curve.*) (v0.75, Research OS R10)
/ ----------------------------------------------------------------------------
/ ARCHITECTURE.md Part II 11.8(b) / 13-R10. The SECOND evidence-layer milestone: a dedicated
/ curve capability that consumes R9's single door (.state.*) and computes the RICH derived
/ curve - roll yield, slope, curvature, contango/backwardation classification - plus the
/ curve-shock operators and immutable daily snapshots. It is the curve source the replay engine
/ (R12), seasonality/carry (R14), PnL/risk attribution (R15) and the first real strategy (R16)
/ all build on; because it reads through R9's door, every curve it produces is point-in-time by
/ construction.
/ -
/ LOW layer: loads AFTER state/ (consumes .state.asof / .state.build) and BEFORE the consumers.
/ Reads state/ (and via it data/) downward only; opens no HDB/snapshots at import. ADDITIVE: it
/ does NOT rewire regime/ or .state.build (regime/ keeps its own slope/classification; rebasing
/ it onto the curve engine is a later deferral). Registered as an R2 `curve capability; carded.
/ -
/ AGREEMENT WITH regime/: the slope and the backwardation/contango/flat classification use the
/ SAME convention + threshold as regime/ (slope = (deferred-front)/front, deferredIdx and the flat
/ threshold sourced from .cfg.regime via .cfg.curve), so the two derivations agree on the same
/ curve - the precondition for eventually rebasing regime/ onto this engine. Proven by a test.
/ -
/ RESERVED-NAME NOTE: parameter is `asOf` NEVER `asof`; `comm` NOT `commodity` (the qSQL
/ self-compare trap). Builtins not shadowed (deltas/avg/dev/var/cor/cov/sum/first/last/differ).
/ ----------------------------------------------------------------------------

.curve.defaultConfig:{[] .cfg.curve};

/ ── derived features (pure: a curve table -> slope/curvature/rollYield/classification) ──

/ The curve is the .data.hdb.curveAt shape: (tenor [years], price, contractYM, expiry), tenor-sorted.
/ slope = (deferred-front)/front (matches regime/); classification by the regime flat threshold;
/ curvature = the butterfly second difference (near - 2*mid + far); rollYield = annualised
/ ln(Ffront/Fdeferred) (positive in backwardation), annualised by the near->deferred tenor gap.
.curve.__features:{[curve;cfg]
    prices:`float$curve`price;
    tenors:`float$curve`tenor;
    n:count prices;
    front:first prices;
    defIdx:(cfg`deferredIdx)&n-1;
    deferred:prices defIdx;
    slope:$[front=0f; 0n; (deferred-front)%front];
    flatThr:cfg`flatThreshold;
    classification:?[slope<neg flatThr;`backwardation;?[slope>flatThr;`contango;`flat]];
    midIdx:n div 2;
    curvature:$[n>=3; (first[prices] + last prices) - 2f * prices midIdx; 0n];   / near - 2*mid + far
    dt:(tenors defIdx)-first tenors;
    rollYield:$[(front>0f) and (deferred>0f) and dt>0f; (log front%deferred)%dt; 0n];
    `slope`curvature`rollYield`classification!(slope;curvature;rollYield;classification)
 };

/ Adjacent calendar spreads from a curve.
.curve.__spreads:{[curve]
    $[1<count curve;
        ([] nearTenor:(-1_curve`tenor); farTenor:(1_curve`tenor);
            nearContract:(-1_curve`contractYM); farContract:(1_curve`contractYM);
            spread:(-1_curve`price)-1_curve`price);
        ([] nearTenor:`float$(); farTenor:`float$(); nearContract:`long$(); farContract:`long$(); spread:`float$())]
 };

/ ── the curve engine (.curve.build) ──────────────────────────────────────────

/ Build the rich curve object for (asOf, comm) THROUGH R9's door, with a liquidity filter
/ (drop non-positive prices + contracts below the volume floor per .cfg.curve). Returns a dict:
/ asOf, commodity, effectiveDate, curve (clean), spreads, features (slope/curvature/rollYield/
/ classification). Deterministic. Reuses .state.build's curve shape (consistent with regime/).
.curve.build:{[asOf;comm]
    st:.state.build[asOf;comm];
    effDate:st`effectiveDate;
    cfg:.cfg.curve;
    / the liquid contracts on the effective date (via the same as-of door), for the tenor filter.
    effRows:select from (.state.asof[asOf;comm])`data where date=effDate;
    liquid:exec contractYM from effRows where volume>=cfg`minVolume;
    curve0:st`curve;
    curve:`tenor xasc select from curve0 where price>0f, contractYM in liquid;
    if[0=count curve; '"curve.build: no liquid tenors for ",(string comm)," as of ",string asOf];
    `asOf`commodity`effectiveDate`curve`spreads`features!(
        asOf; comm; effDate; curve; .curve.__spreads curve; .curve.__features[curve;cfg])
 };

/ ── curve-shock operators (pure: (curve;size) -> shockedCurve) ───────────────

/ Parallel: add a constant to every tenor.
.curve.shock.parallel:{[curve;size] update price:price+size from curve};
/ Slope: add a tenor-proportional amount, pivoting on the front (front unchanged, far moves most).
.curve.shock.slope:{[curve;size] update price:price+size*tenor-first tenor from curve};
/ Butterfly: add a curvature-shaped amount (a 4x(1-x) bow over the tenor span - belly moves, ends fixed).
.curve.shock.butterfly:{[curve;size]
    t:`float$curve`tenor;
    x:$[(max t)>min t; (t-min t)%(max t)-min t; 0f*t];
    bow:4f*x*1f-x;
    update price:price+size*bow from curve};

/ ── immutable daily snapshots (splayed, UNPARTITIONED - matching the HDB decision) ──

.curve.__emptySnap:{[]
    ([] commodity:`symbol$(); date:`date$(); asOf:`date$();
        slope:`float$(); curvature:`float$(); rollYield:`float$(); classification:`symbol$(); nTenors:`long$())
 };
.curve.__initSnap:{[] if[not `snapshots in key `.curve; .curve.snapshots:.curve.__emptySnap[]];};

/ A snapshot row from a built curve object.
.curve.__snapRow:{[c]
    f:c`features;
    `commodity`date`asOf`slope`curvature`rollYield`classification`nTenors!(
        c`commodity; c`effectiveDate; c`asOf; f`slope; f`curvature; f`rollYield; f`classification; count c`curve)
 };

/ Immutable write: a snapshot for (commodity, date) is written ONCE and NEVER overwritten in
/ place. Recomputing must produce the identical curve (deterministic), so a re-snapshot is a
/ no-op-if-matching; if it would DIFFER it errors (catching non-determinism). Compares on the
/ feature fields (asOf may legitimately differ for a later as-of of the same observation date,
/ so the canonical key + match is on (commodity, date) + the curve features).
.curve.__snapWrite:{[row]
    .curve.__initSnap[];
    cmp:`slope`curvature`rollYield`classification`nTenors;
    existing:select from .curve.snapshots where commodity=row`commodity, date=row`date;
    if[count existing;
        if[not (cmp#first existing)~cmp#row;
            '"curve.snapshot: a DIFFERING snapshot already exists for ",(string row`commodity),"/",(string row`date)," (non-deterministic curve)"];
        :row];                              / no-op if matching (immutable)
    .curve.snapshots:.curve.snapshots upsert row;
    row
 };

/ Compute via .curve.build and persist immutably to the in-memory snapshots table.
.curve.snapshot:{[asOf;comm] .curve.__snapWrite .curve.__snapRow .curve.build[asOf;comm]};
.curve.snapshotAt:{[date;comm] .curve.__initSnap[]; select from .curve.snapshots where commodity=comm, date=date};
.curve.snapshotHistory:{[comm] .curve.__initSnap[]; `date xasc select from .curve.snapshots where commodity=comm};

/ Persist / load the snapshots splay (splayed, UNPARTITIONED, p# on commodity, .Q.en + set -
/ the SAME storage decision as the HDB; NOT date-partitioned). gitignored; tests synthetic only.
.curve.snapshotPersist:{[hdbPath]
    .curve.__initSnap[];
    t:.Q.en[hsym `$hdbPath; update `p#commodity from `commodity`date xasc .curve.snapshots];
    (hsym `$hdbPath,"/curveSnapshots/") set t;
    `rows!enlist count t
 };
.curve.snapshotOpen:{[hdbPath]
    sp:hsym `$hdbPath,"/curveSnapshots/";
    if[0=count key sp; '"curve.snapshotOpen: no curveSnapshots at ",hdbPath];
    .curve.snapshots:get sp;
    hdbPath
 };

/ ── register as an R2 `curve capability ──────────────────────────────────────

.curve.contract:`version`requiredIn`requiredOut!(
    1;
    `asOf`commodity!`date`symbol;
    `asOf`commodity`curve`spreads`features!`date`symbol`table`table`dict);
.registry.new[`curve; .curve.contract];
.curve.register:{[name;fn;manifest] .registry.register[`curve;name;fn;manifest]};
.curve.get:{[name] .registry.get[`curve;name]};
.curve.list:{[] .registry.list `curve};
.curve.conforms:{[name] .registry.conforms[`curve;name]};

.curve.register[`curveEngine; .curve.build;
    `contractVersion`description`in`out!(
        1;
        "the curve engine: clean as-of curve + spreads + derived features (roll yield/slope/curvature/classification)";
        `asOf`commodity!`date`symbol;
        `asOf`commodity`curve`spreads`features!`date`symbol`table`table`dict)];

.curve.__initSnap[];
-1 "curve.q loaded - .curve.* curve engine + shocks + immutable snapshots ready (not opened)";
