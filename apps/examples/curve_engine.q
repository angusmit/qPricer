\l core/init.q
/ ============================================================================
/ curve_engine.q - real-data example (NOT a test), DATA-CONDITIONAL. Research OS R10.
/ Build the rich curve for real CRUDE through R9's as-of door: print slope / curvature /
/ roll yield / classification; show a backwardated date vs a contango date (the classification
/ flips); apply the three shock operators; snapshot a few dates and read them back (round-trip
/ + immutability). Real crude. End exit 0;.
/ ============================================================================
hdbPath:.cfg.paths`hdb;
useHdb:0<count @[{[p] key hsym `$p,"/sym"};hdbPath;{[e] ()}];
if[not useHdb;
    -1 "curve_engine SKIPPED - no HDB at ",hdbPath," (run scripts/ingest_hdb.q first).";
    exit 0];
.data.hdb.open hdbPath;

dts:.data.hdb.dates `CRUDE;
/ sample ~24 dates evenly, build each, and find one backwardation + one contango date.
sample:dts where 0=(til count dts) mod (1|(count dts) div 24);
built:{[d] @[{[d] c:.curve.build[d;`CRUDE]; `date`class`slope`curvature`roll!(d;(c`features)`classification;(c`features)`slope;(c`features)`curvature;(c`features)`rollYield)};d;{[e] ()}]} each sample;
built:built where 0<count each built;
backDate:first exec date from built where class=`backwardation;
contDate:first exec date from built where class=`contango;

showCurve:{[d]
    c:.curve.build[d;`CRUDE]; f:c`features;
    -1 "  as of ",(string d)," (effective ",(string c`effectiveDate),"): ",(string f`classification),
        "  slope=",(string f`slope),"  curvature=",(string f`curvature),"  rollYield=",string f`rollYield;
    -1 "    front 3 tenors: ",-3!3 sublist c`curve};

-1 "RICH CURVE (CRUDE) - the classification flips between regimes:";
$[null backDate; -1 "  (no backwardated date in the sample)"; showCurve backDate];
$[null contDate; -1 "  (no contango date in the sample)"; showCurve contDate];
-1 "";

/ the three shock operators on a real curve.
refDate:$[not null contDate; contDate; first exec date from built];
base:(.curve.build[refDate;`CRUDE])`curve;
-1 "Curve-shock operators on the ",(string refDate)," curve (front 4 prices):";
-1 "  base       : ",-3!4 sublist base`price;
-1 "  parallel+1 : ",-3!4 sublist (.curve.shock.parallel[base;1f])`price;
-1 "  slope+1    : ",-3!4 sublist (.curve.shock.slope[base;1f])`price;
-1 "  butterfly+1: ",-3!4 sublist (.curve.shock.butterfly[base;1f])`price;
-1 "";

/ snapshot a few dates, read them back, show immutability.
.curve.snapshots:.curve.__emptySnap[];
snapDates:5 sublist exec date from built;
.curve.snapshot[;`CRUDE] each snapDates;
-1 "Immutable daily snapshots (",(string count snapDates)," dates):";
show .curve.snapshotHistory[`CRUDE];
.curve.snapshot[first snapDates;`CRUDE];   / re-snapshot = no-op
-1 "after re-snapshotting the first date, history still has ",(string count .curve.snapshotHistory[`CRUDE])," rows (immutable)";
-1 "";
-1 "(R10: the curve source the evidence layer composes - point-in-time by construction (built through";
-1 " R9's door), agreeing with regime/'s classification convention, with shock operators for R15 risk";
-1 " attribution and immutable snapshots that make a replayed curve reproducible rather than recomputed.)";
exit 0;
