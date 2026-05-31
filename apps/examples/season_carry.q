\l core/init.q
/ ============================================================================
/ season_carry.q - real-data example (NOT a test), DATA-CONDITIONAL.
/ Research OS R14 (ARCHITECTURE.md 11.8(f)): the evidence-layer FEATURE capabilities.
/ SEASONALITY - the CAUSAL same-calendar-month z-score of the real CRUDE front-deferred spread on a
/ chosen date vs its same-month history UP TO that date; contrast with the FULL-SAMPLE z (which would
/ use FUTURE same-months - the look-ahead a seasonal backtest must avoid). CARRY - the implied carry /
/ convenience yield (flagged assumption-dependent) / carry signal on a backwardated date vs a contango
/ date (the carry signal flips sign). Skips gracefully if the HDB is absent.
/ ============================================================================
hdbPath:.cfg.paths`hdb;
useHdb:0<count @[{[p] key hsym `$p,"/sym"};hdbPath;{[e] ()}];
if[not useHdb;
    -1 "season_carry SKIPPED - no HDB at ",hdbPath," (run scripts/ingest_hdb.q first).";
    exit 0];
.data.hdb.open hdbPath;

comm:`CRUDE;
dts:.data.hdb.dates comm;
asOf:last dts;

/ (1) SEASONALITY - causal same-month z vs a full-sample (look-ahead) z.
-1 "(1) Seasonality - causal same-calendar-month z of the CRUDE front-deferred spread as of ",(string asOf),":";
f:.season.features[asOf;comm];
-1 "    spread=",(string f`spread),"  sameMonthZ(causal)=",(string f`sameMonthZ),
   "  nSameMonth=",(string f`nSameMonth),"  lowConfidence=",string f`lowConfidence;
-1 "    sameContractMonthZ=",(string f`sameContractMonthZ),"  seasonalFactor=",(string f`seasonalFactor),
   "  deseasonalised=",string f`deseasonalised;
/ the LOOK-AHEAD contrast: a full-sample same-month z would use the WHOLE series (incl. dates AFTER asOf).
ser:.season.__series[asOf;comm;.cfg.season`deferredIdx];        / causal series (date<=asOf)
serFull:.season.__series[last dts;comm;.cfg.season`deferredIdx]; / the full series (here last dts == asOf, illustrative)
cm:(`month$ser`date) mod 12; curCm:(`month$asOf) mod 12;
causalSm:(ser`spread) where cm=curCm;
-1 "    -> causal same-month uses ",(string count causalSm)," observations, ALL dated <= ",(string asOf);
-1 "       (a full-sample seasonal stat would standardise this past spread by FUTURE same-months - look-ahead).";
-1 "";

/ (2) CARRY - on a backwardated date vs a contango date (the carry signal flips sign).
-1 "(2) Carry - implied carry / convenience yield (assumption-dependent) / carry signal:";
labs:.regime.series[comm;dts];
cs:labs`curveState;                                        / curveState per date (aligned to dts)
cleanIdx:where dts<(last dts)-35;                          / avoid the trailing data-edge (sparse curve)
bwAll:dts cleanIdx where `backwardation=cs cleanIdx;
coAll:dts cleanIdx where `contango=cs cleanIdx;
/ pick the MEDIAN date of each regime (the heart of the regime - a clean multi-contract curve, away
/ from the sparse data-edge tail where the deferred leg jumps to a far contract).
bwDate:$[count bwAll; bwAll (count bwAll) div 2; 0Nd];
coDate:$[count coAll; coAll (count coAll) div 2; 0Nd];
showCarry:{[comm;d;tag] c:.carry.features[d;comm];
    -1 "  ",tag," ",(string d),": class=",(string c`classification),
       "  impliedCarry=",(string c`impliedCarry),"  convYield=",(string c`convenienceYield),
       "  carrySignal=",string c`carrySignal};
if[not null bwDate; showCarry[comm;bwDate;"backwardation"]];
if[not null coDate; showCarry[comm;coDate;"contango     "]];
-1 "    (impliedCarry + the carry signal are POSITIVE in backwardation, NEGATIVE in contango - they flip;";
-1 "     convenience yield + cash-and-carry fair value depend on the ASSUMED r+storage from .cfg.carry.)";

exit 0;
