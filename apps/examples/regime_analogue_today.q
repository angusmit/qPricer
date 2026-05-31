\l core/init.q
/ ============================================================================
/ regime_analogue_today.q - real-data example (NOT a test), DATA-CONDITIONAL.
/ Research OS R4 (ARCHITECTURE.md 11.3): take the LATEST crude date's regime
/ fingerprint and ask the analogue engine "this state most resembles WHAT, and what
/ killed strategies there?" - the digested-history payoff, explicit and queryable.
/ -
/ Needs the HDB + scripts/build_regimes.q + scripts/build_regime_library.q to have run.
/ Skips gracefully otherwise. End exit 0;.
/ ============================================================================
hdbPath:.cfg.paths`hdb;
useHdb:0<count @[{[p] key hsym `$p,"/sym"};hdbPath;{[e] ()}];
if[not useHdb;
    -1 "regime_analogue_today SKIPPED - no HDB at ",hdbPath," (run scripts/ingest_hdb.q first).";
    exit 0];
.data.hdb.open hdbPath;
if[0=count key hsym `$hdbPath,"/regimeEpisodes/";
    -1 "regime_analogue_today SKIPPED - no regimeEpisodes (run scripts/build_regimes.q then scripts/build_regime_library.q).";
    exit 0];
.regime.library.open hdbPath;

today:last .data.hdb.dates `CRUDE;
state:.regime.label[`CRUDE;today];
-1 "Crude regime fingerprint on ",(string today),":";
-1 "  curveState=",(string state`curveState),"  volState=",(string state`volState),"  liqState=",(string state`liqState),
   "  rollPhase=",(string state`rollPhase),"  seasonPhase=",string state`seasonPhase;
-1 "  slopePct=",(string state`slopePct),"  volPct=",(string state`volPct),"  volumePct=",string state`volumePct;
-1 "";

an:.regime.analogue.forDate[`CRUDE;today;3];
-1 "Today's crude state most resembles (nearest historical episodes):";
{[r] -1 "";
    -1 "  ",(string r`label)," (",(string r`dateFrom),"..",(string r`dateTo),")  distance=",string r`distance;
    -1 "    drivers: ",string r`driversKey;
    -1 "    RISK MEMORY: ",r`riskMemory} each an;
-1 "";
-1 "(This is what an experienced desk does in their head - 'this resembles 2020 / 2022' and 'here is what";
-1 " blew up then' - now an explicit, queryable record rather than a recollection. Pair with the R3/R3b";
-1 " gates: gov MAY consult .regime.analogue (a legal downward call) to surface these failure modes when";
-1 " gating a hypothesis.)";
exit 0;
