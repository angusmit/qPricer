\l core/init.q
/ ============================================================================
/ relative_value_template.q - real-data example (NOT a test), DATA-CONDITIONAL.
/ Research OS R6: a NEW research SHAPE (relative-value) through the SAME governance.
/ Instantiate the relativeValue template on a concrete crude CALENDAR spread (front vs
/ the next tenor), print the spread + mean-reversion signal, run its TEMPLATE-SPECIFIC
/ spread-stationarity gate, then - only the universal part - flow the PnL through
/ .gov.runFull (cost -> deflation -> walk-forward -> sealed holdout). Honest verdict; not
/ tuned to pass. Skips if no HDB. exit 0;.
/ ============================================================================
hdbPath:.cfg.paths`hdb;
useHdb:0<count @[{[p] key hsym `$p,"/sym"};hdbPath;{[e] ()}];
if[not useHdb;
    -1 "relative_value_template SKIPPED - no HDB at ",hdbPath," (run scripts/ingest_hdb.q first).";
    exit 0];
.data.hdb.open hdbPath;

clHist:.data.hdb.curveHistory[`CRUDE;.data.hdb.dates `CRUDE];
inputs:`commodity`curveHistory!(`CRUDE;clHist);

/ --- the template-specific gate runs FIRST (a gate directional does not have) ---
full:.template.run[`relativeValue; inputs];
stat:full[`validation;`stationarity];
-1 "relativeValue on a CRUDE calendar spread (front vs next tenor), ",(string count full`pnl)," days";
-1 "spread-stationarity gate (TEMPLATE-SPECIFIC): AR(1)=",(string stat`ar1),"  half-life=",(string stat`halfLife)," days  stationary=",string stat`stationary;
-1 "spread PnL: total=",(string sum (full`pnl)`pnl),"  (edge source: structural / spread reversion)";
-1 "";

if[not stat`stationary;
    -1 "The spread FAILS the stationarity gate -> the RV thesis is rejected BEFORE the universal";
    -1 "gov gates even run (a non-stationary / random-walk spread has no mean-reversion edge).";
    -1 "This is the template carrying its OWN gate that the directional template does not have."];

if[stat`stationary;
    -1 "The spread is stationary -> proceed to the UNIVERSAL gov gates via .gov.runFull.";
    runner:{[inputs;from;to] (.template.rv.run[inputs,`dateFrom`dateTo!(from;to)])`pnl}[inputs];
    .gov.hypoTbl:.gov.__emptyHypotheses[]; .gov.trialTbl:.gov.__emptyTrials[];
    .gov.register `hypoId`thesis`edgeSource`instruments`claimedRegimes`status!(
        `rv_crude_cal;"crude calendar spread mean-reverts (structural)";`mispricing;`CRUDE;`contango`backwardation;`research);
    verdicts:.gov.runFull[`rv_crude_cal;runner;`curveState];
    -1 "";
    -1 "Universal gov cascade by curve state (the SAME gates the directional shapes face):";
    show select bucket,nObs,netSharpe,dsr,failedGate,verdict,tradeable from verdicts];

-1 "";
-1 "(R6: a genuinely different SHAPE - composed from capabilities, carrying its OWN stationarity";
-1 " gate - judged by the same honest gov gates. The spine now holds more than directional bets.)";
exit 0;
