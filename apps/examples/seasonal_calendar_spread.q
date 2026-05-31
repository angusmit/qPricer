\l core/init.q
/ ============================================================================
/ seasonal_calendar_spread.q - real-data example (NOT a test), DATA-CONDITIONAL. THE R16 CAPSTONE.
/ Run seasonally-adjusted crude calendar-spread mean-reversion END-TO-END through the WHOLE foundation
/ (R9 door -> R10 curve -> R11 roll -> R12 replay -> R13 evidence audit -> R14 seasonality -> R15
/ attribution) via .workflow.runReplay, and print the HONEST verdict + the PnL attribution + the
/ human-escalation packet. The parameters are PRE-REGISTERED (.cfg.strategy.seasonalCalSpread) and NOT
/ tuned to pass: a "research, not tradeable" verdict is the system WORKING, not a failure. Skips if no HDB.
/ ============================================================================
hdbPath:.cfg.paths`hdb;
useHdb:0<count @[{[p] key hsym `$p,"/sym"};hdbPath;{[e] ()}];
if[not useHdb;
    -1 "seasonal_calendar_spread SKIPPED - no HDB at ",hdbPath," (run scripts/ingest_hdb.q first).";
    exit 0];
.data.hdb.open hdbPath;

comm:`CRUDE; cfg:.cfg.strategy.seasonalCalSpread;
/ scope the audited run record to train+validate (the holdout stays sealed; gov reads it one-shot).
tv:.gov.zone.range[comm;`trainValidate];
res:.template.scs.run `commodity`dateFrom`dateTo!(comm; tv 0; tv 1);
run:(res`meta)`runRecord;
-1 "seasonalCalendarSpread on ",(string comm)," (train+validate ",(string tv 0)," .. ",(string tv 1),"):";
-1 "  stationarity: ",(res[`validation;`stationarity])`detail;
-1 "  run record: ",(string count run`steps)," steps, replay totalPnl=",string (run`meta)`totalPnl;
-1 "";

/ run the END-TO-END bounded replay-mode loop: carded -> evidence audit -> gates -> skeptic -> attribution.
/ NOTE: the gov hypothesis declares edgeSource `mispricing (gov's validEdgeSources vocabulary for the
/ structural relative-value / dislocation edge); the model card uses the cards-layer term `structural -
/ the same edge, two enums (gov: riskPremium/mispricing/info; cards: riskPremium/structural/informational).
prop:`hypoId`thesis`edgeSource`instruments`claimedRegimes`capName`runRecord`runner!(
    `crude_seasonalCalSpread;
    "the seasonally-adjusted crude M1-M3 calendar spread reverts to its seasonal norm (relative-value/mispricing)";
    `mispricing; enlist comm; `backwardation`contango`flat; `seasonalCalSpread;
    run; .template.scs.runner[comm;cfg]);
pkt:.workflow.runReplay prop;

-1 "(1) The HONEST verdict (per regime bucket; gates: cost -> deflated Sharpe -> walk-forward -> sealed holdout):";
$[98h=type pkt`verdicts;
    show select bucket,nObs,verdict,failedGate,tradeable,netSharpe,dsr from pkt`verdicts;
    -1 "    ",-3!pkt`verdicts];
-1 "  anyTradeable = ",(string pkt`anyTradeable),"   (evidence audit passed = ",(string pkt`evidencePassed),")";
-1 "  regime risk-memory: ",$[0<count pkt`riskMemory; pkt`riskMemory; "(none)"];
-1 "";

-1 "(2) PnL ATTRIBUTION (R15: where did the return come from?):";
a:pkt`attribution;
-1 "    level=",(string a`level),"  slope=",(string a`slope),"  curvature=",(string a`curvature),"  carry=",(string a`carry),"  residual=",string a`residual;
-1 "    total=",(string a`total),"  RESIDUAL FRACTION=",string a`residualFraction;
-1 "    (note: a single-contract attribution captures the FRONT leg; the deferred leg of a 2-leg spread";
-1 "     appears in residual - a known limitation, not necessarily noise. A genuinely large residual on a";
-1 "     single-leg run would be the red flag that the edge is not a recognisable curve/carry source.)";
-1 "";

-1 "(3) The human-escalation packet:";
-1 "    decision=",(string pkt`decision),"  humanSignOffRequired=",(string pkt`humanSignOffRequired),"  stage=",string pkt`stage;
-1 "    summary: ",pkt`summary;
-1 "";
-1 "This is the capstone: the verdict is HONEST - judged on realistic replay evidence (R12) that passed";
-1 "the evidence audit (R13), deflated against the trial count (R3), tested out-of-sample (walk-forward +";
-1 "the sealed one-shot holdout), stressed by the regime skeptic, and attributed (R15). A not-tradeable";
-1 "verdict is the architecture doing its job: telling the truth instead of a flattering backtest.";

exit 0;
