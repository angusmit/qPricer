\l core/init.q
/ ============================================================================
/ pnl_attribution.q - real-data example (NOT a test), DATA-CONDITIONAL.
/ Research OS R15 (ARCHITECTURE.md 11.8(g)): turn "explain why this makes money" into a QUANTITATIVE
/ attribution. Take a real CRUDE replay run (R12 momentum), decompose its realized PnL into
/ level / slope / curvature / carry / residual (R10's shock basis), show it RECONCILES to the realized
/ total, and print the RESIDUAL FRACTION (a large residual = the edge is NOT a recognisable curve/carry
/ source - the quantitative "name which of the three edges, or it's noise"). Then the bucketed curve
/ risk: the per-tenor delta + roll-down + calendar-spread exposure. Skips if no HDB.
/ ============================================================================
hdbPath:.cfg.paths`hdb;
useHdb:0<count @[{[p] key hsym `$p,"/sym"};hdbPath;{[e] ()}];
if[not useHdb;
    -1 "pnl_attribution SKIPPED - no HDB at ",hdbPath," (run scripts/ingest_hdb.q first).";
    exit 0];
.data.hdb.open hdbPath;

comm:`CRUDE; strat:`timeSeriesMomentum;
dts:.data.hdb.dates comm;
/ a ~120-date window inside the clean curve range (avoid the sparse data-edge tail).
toD:dts[(count dts)-40 & count dts];
fromD:dts[(count dts)-160 & count dts];
run:.backtest.replay.run[strat;comm;fromD;toD;()!()];
-1 "Replay run: ",(string strat)," on ",(string comm),", ",(string fromD)," .. ",(string toD),
   "  (realized totalPnl=",(string (run`meta)`totalPnl),")";
-1 "";

/ (1) PnL decomposition.
a:.attribution.pnl run;
-1 "(1) PnL attribution (level / slope / curvature / carry / residual):";
-1 "    level     = ",string a`level;
-1 "    slope     = ",string a`slope;
-1 "    curvature = ",string a`curvature;
-1 "    carry     = ",string a`carry;
-1 "    residual  = ",string a`residual;
-1 "    --------------------------------";
-1 "    total     = ",(string a`total),"   (reconciled=",(string a`reconciled),")";
-1 "    RESIDUAL FRACTION = ",string a`residualFraction;
-1 "    -> a large residual fraction means the PnL is NOT explained by a recognisable curve/carry";
-1 "       source - the quantitative red flag that the edge may be noise / a missing factor.";
-1 "";

/ (2) bucketed curve risk.
rk:.attribution.risk run;
-1 "(2) Bucketed curve risk on the end-of-run position (",(string rk`endContract)," @ tenor ",(string rk`tenor),"yr):";
-1 "    bucketDelta (per tenor bucket ",(" " sv string rk`buckets),"):";
show rk`bucketDelta;
-1 "    rollDownExposure       = ",string rk`rollDownExposure;
-1 "    calendarSpreadExposure = ",string rk`calendarSpreadExposure;

exit 0;
