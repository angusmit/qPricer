\l core/init.q
/ ============================================================================
/ market_state_asof.q - real-data example (NOT a test), DATA-CONDITIONAL. Research OS R9.
/ The single door to history, demonstrated: build the point-in-time Market State for a real
/ CRUDE as-of date (print the as-of curve / calendar spreads / tradable universe / provenance),
/ then show that an EARLIER as-of date is strictly a time-prefix - it excludes the later
/ contracts and later prices. Real crude. End exit 0;.
/ ============================================================================
hdbPath:.cfg.paths`hdb;
useHdb:0<count @[{[p] key hsym `$p,"/sym"};hdbPath;{[e] ()}];
if[not useHdb;
    -1 "market_state_asof SKIPPED - no HDB at ",hdbPath," (run scripts/ingest_hdb.q first).";
    exit 0];
.data.hdb.open hdbPath;

dts:.data.hdb.dates `CRUDE;
/ use a late as-of ~60 trading days before the data ends (NOT the final date - there the HDB
/ expiry, derived as each contract's last-observed date, equals the as-of, so no contract reads
/ as live; that boundary is honest given the data but a poor illustration).
asOfLate:dts (count[dts]-60)&(count dts)-1;
asOfEarly:dts (count dts) div 3;       / an as-of date ~1/3 through history

stLate:.state.build[asOfLate;`CRUDE];
stEarly:.state.build[asOfEarly;`CRUDE];

-1 "MARKET STATE for CRUDE as of ",(string asOfLate)," (effective date ",(string stLate`effectiveDate),"):";
-1 "  as-of curve (the same shape regime/ sees):";
show stLate`curve;
-1 "  calendar spreads (adjacent maturities):";
show stLate`spreads;
-1 "  tradable universe (contracts not expired as of asOf): ",(string count stLate`universe)," contracts ",-3!stLate`universe;
-1 "  refs: ",-3!stLate`refs;
-1 "  provenance: as-of ",(string (stLate`provenance)`asOf),", ",(string (stLate`provenance)`nRows)," rows, dates ",(string (stLate`provenance)`dateFrom)," .. ",string (stLate`provenance)`dateTo;
-1 "";

-1 "AS-OF DISCIPLINE - the earlier as-of date is strictly a time-prefix:";
-1 "  early as-of ",(string asOfEarly),": effDate=",(string stEarly`effectiveDate),", ",(string (stEarly`provenance)`nRows)," rows, ",(string count stEarly`universe)," tradable contracts";
-1 "  late  as-of ",(string asOfLate),": effDate=",(string stLate`effectiveDate),", ",(string (stLate`provenance)`nRows)," rows, ",(string count stLate`universe)," tradable contracts";
-1 "  early effDate strictly before late: ",string (stEarly`effectiveDate) < stLate`effectiveDate;
-1 "  early sees no row after its as-of: ",string .state.invariant.asofRespected[asOfEarly; (.state.asof[asOfEarly;`CRUDE])`data];
-1 "  later contracts unknown early: the early universe is a subset in time (max contract ",(string max stEarly`universe)," <= late max ",(string max stLate`universe),"): ",string (max stEarly`universe)<=max stLate`universe;
-1 "";
-1 "(R9: ONE as-of door, ONE Market State object. Point-in-time is now enforceable by construction -";
-1 " the accessor returns only what is knowable as of asOf and logs provenance; R13's evidence audit";
-1 " will VERIFY that provenance. Every later foundation milestone (curve R10, roll R11, replay R12,";
-1 " season/carry R14, the first real strategy R16) consumes this state instead of raw prices.)";
exit 0;
