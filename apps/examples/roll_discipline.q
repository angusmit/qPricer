\l core/init.q
/ ============================================================================
/ roll_discipline.q - real-data example (NOT a test), DATA-CONDITIONAL.
/ Research OS R11 (ARCHITECTURE.md 11.8(c)): roll discipline on real CRUDE. Shows (1) which
/ contract is ACTIVE as of a date, decided FROM AS-OF DATA ONLY (through R9's door); (2) the
/ roll events over a window (front -> next, with the reason); (3) the analytics-only continuous
/ back-adjusted series vs the ACTUAL contract you would trade; (4) the as-of-only property - a
/ past decision is fixed by the data knowable then, not by anything later. Skips if no HDB.
/ THE PRINCIPLE: trade the actual contract the roll map names; never trade/signal off continuous.
/ ============================================================================
hdbPath:.cfg.paths`hdb;
useHdb:0<count @[{[p] key hsym `$p,"/sym"};hdbPath;{[e] ()}];
if[not useHdb;
    -1 "roll_discipline SKIPPED - no HDB at ",hdbPath," (run scripts/ingest_hdb.q first).";
    exit 0];
.data.hdb.open hdbPath;

comm:`CRUDE;
dts:.data.hdb.dates comm;
/ pick the latest date that actually has a live (un-expired) deferred contract - the clean curve
/ range (cheap, from `futures`; this only CHOOSES where to point the demo, not a roll decision).
me:select maxExp:max expiry by date from futures where commodity=comm;
liveDates:asc exec date from me where maxExp>date;
if[0=count liveDates; -1 "roll_discipline SKIPPED - no live deferred contracts in the CRUDE data."; exit 0];
-1 "Roll discipline on ",(string comm)," - rule: ",string (.roll.__rule comm)`type;
-1 "";

/ (1) the ACTIVE contract as of a recent date, from as-of data only.
asOf:last liveDates;
act:.roll.active[asOf;comm];
-1 "Active contract as of ",(string asOf),": ",(string act`active),
   "  (reason=",(string act`reason),", effectiveDate=",(string act`effectiveDate),")";
-1 "  weights: ",(", " sv {[k;v] (string k),"=",string v}'[key act`weights; value act`weights]);
-1 "";

/ (2) roll events over the last ~year of the clean curve range.
fromDate:liveDates (count liveDates)-260 & count liveDates;
ev:.roll.events[comm;fromDate;asOf];
-1 "Roll events ",(string fromDate)," -> ",(string asOf)," (",(string count ev)," rolls):";
show ev;
-1 "";

/ (3) the continuous back-adjusted series (ANALYTICS-ONLY) vs the ACTUAL contract price.
cont:.roll.continuous[comm];
tail:-8#cont;
-1 "Continuous (back-adjusted, analytics-only) vs the actual active-contract front price:";
show select date, activeContract:contract, actualFront:rawFront, continuousAdj:adjusted, rolled from tail;
-1 "  NOTE: 'continuousAdj' rewrites past levels at each roll (non-point-in-time) - chart/vol ONLY.";
-1 "  'actualFront' is the price of the contract .roll.active names - the tradable, point-in-time path.";
-1 "";

/ (4) the as-of-only property: the decision at a PAST date is the same whether or not later data exists.
pastDate:liveDates (count liveDates)-30 & (count liveDates)-1;
full:(.roll.active[pastDate;comm])`active;
-1 "As-of-only property: active as of ",(string pastDate)," = ",(string full),
   " - decided from data with date<=",(string pastDate)," ONLY (no future bar can change it).";

exit 0;
