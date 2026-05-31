\l core/init.q
/ ============================================================================
/ replay_engine.q - real-data example (NOT a test), DATA-CONDITIONAL.
/ Research OS R12 (ARCHITECTURE.md 11.8(d)): THE KEYSTONE. Replay a real COMMODITY futures strategy
/ (cross-commodity momentum on CRUDE) through the event-driven, as-of-verified replay loop and show:
/ (1) FAITHFULNESS - a frictionless replay AGREES with the research-mode result to a tight tolerance
/     (the look-ahead detector: it reproduces the vectorised total when there is no look-ahead);
/ (2) REALISM - turn ON the participation cap + roll-window penalty and watch the net PnL DEGRADE
/     (gross vs net, turnover, cost drag - the realistic evidence the gates will judge);
/ (3) the auditable RUN RECORD - the book, a few fills, the roll events, the per-step PnL components
/     (price / cost / financing) + the as-of provenance (R13's evidence audit bites on exactly this).
/ Skips gracefully if the HDB is absent.
/ ============================================================================
hdbPath:.cfg.paths`hdb;
useHdb:0<count @[{[p] key hsym `$p,"/sym"};hdbPath;{[e] ()}];
if[not useHdb;
    -1 "replay_engine SKIPPED - no HDB at ",hdbPath," (run scripts/ingest_hdb.q first).";
    exit 0];
.data.hdb.open hdbPath;

strat:`timeSeriesMomentum; comm:`CRUDE;

/ (1) FAITHFULNESS - frictionless replay vs research mode (the same strategy, untouched research path).
-1 "(1) Faithfulness - frictionless replay vs research-mode PnL (full CRUDE history):";
f:.backtest.replay.faithfulness[strat;comm;1e-6];
-1 "    research totalPnl = ",string f`researchTotal;
-1 "    replay   totalPnl = ",string f`replayTotal;
-1 "    relative diff = ",(string f`relDiff),"   pass=",(string f`pass),"  tol=",(string f`tol),"  (data-edge dates excluded: ",(string f`nEdge),")";
-1 "    -> a frictionless replay reproduces the research total on every date where the roll map agrees;";
-1 "       a divergence there would mean look-ahead or a roll/fill bug (the detector). The few trailing";
-1 "       data-edge dates (research holds a delisted contract) are reported, not hidden.";
-1 "";

/ (2) REALISM - replay a recent window frictionless, then with the participation cap + roll penalty.
dts:.data.hdb.dates comm;
fromD:dts[(count dts)-252 & count dts];
rrF:.backtest.replay.run[strat;comm;fromD;last dts;()!()];
realism:enlist[`replay]!enlist `participationRate`rollPenaltyBps`remainderPolicy!(0.03;8f;`carry);
rrR:.backtest.replay.run[strat;comm;fromD;last dts;realism];
-1 "(2) Realism on the last ~252 dates - frictionless (gross) vs participation cap + roll penalty (net):";
-1 "    frictionless  totalPnl = ",(string (rrF`meta)`totalPnl),"   totalCost = ",(string (rrF`meta)`totalCost);
-1 "    with realism  totalPnl = ",(string (rrR`meta)`totalPnl),"   totalCost = ",(string (rrR`meta)`totalCost);
-1 "    cost drag (gross-net)  = ",(string ((rrF`meta)`totalPnl)-(rrR`meta)`totalPnl),"   (",(string (rrR`meta)`nRolls)," rolls)";
-1 "";

/ (3) the auditable RUN RECORD (from the realism run).
-1 "(3) The auditable run record (R13's evidence audit consumes this):";
-1 "  roll events:"; show rrR`rollEvents;
-1 "  last 6 steps (book + fills + per-step PnL components):";
show -6#select stepDate,activeContract,position,filledQty,refPrice,proportionalCost,slippageCost,rollCost,positionPnl,stepPnl from rrR`steps;
-1 "  provenance (per-step as-of: provDateTo<=asOf always -> no look-ahead):";
show -4#select stepDate,asOf,provNRows,provDateTo from rrR`provenance;

exit 0;
