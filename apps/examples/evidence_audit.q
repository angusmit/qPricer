\l core/init.q
/ ============================================================================
/ evidence_audit.q - real-data example (NOT a test), DATA-CONDITIONAL.
/ Research OS R13 (ARCHITECTURE.md 11.8(e)): the EVIDENCE AUDIT closes the loop the door opened.
/ (1) take a CLEAN real CRUDE replay run (R12) kept inside train+validate -> .evidence.audit -> a PASS
/ report (every check green); (2) deliberately CONTAMINATE it (plant a look-ahead access + a PnL that
/ does not tie) -> .evidence.audit -> a FAIL report (the specific reason); (3) .evidence.gatedRun
/ REJECTS the contaminated run (evidenceFailed, the gates NEVER run) while the clean run is PASSED
/ through to the gates. Skips gracefully if the HDB is absent.
/ ============================================================================
hdbPath:.cfg.paths`hdb;
useHdb:0<count @[{[p] key hsym `$p,"/sym"};hdbPath;{[e] ()}];
if[not useHdb;
    -1 "evidence_audit SKIPPED - no HDB at ",hdbPath," (run scripts/ingest_hdb.q first).";
    exit 0];
.data.hdb.open hdbPath;

comm:`CRUDE; strat:`timeSeriesMomentum;
dts:.data.hdb.dates comm;
boundary:.evidence.__trainValidateBoundary comm;            / the train+validate zone boundary
bi:dts?boundary;
fromD:dts 0|bi-200;                                          / a ~200-date window ending in trainValidate
cleanRun:.backtest.replay.run[strat;comm;fromD;boundary;()!()];

-1 "(1) A clean CRUDE replay run inside train+validate (",(string fromD)," .. ",(string boundary),"):";
audClean:.evidence.audit cleanRun;
show audClean`checks;
-1 "    overall pass = ",(string audClean`pass),"   (",(audClean`reason),")";
-1 "";

/ (2) contaminate: plant a look-ahead access on one step + a PnL that does not reconcile.
s:cleanRun`steps;
pd:s`provDateTo; pd[5]:pd[5]+1;                              / step 5 "saw" data dated after its asOf
sp:s`stepPnl; sp[10]:sp[10]+1.0;                            / step 10's PnL no longer ties
badRun:@[cleanRun;`steps;:;@[@[s;`provDateTo;:;pd];`stepPnl;:;sp]];
-1 "(2) The SAME run with a planted look-ahead + a non-tying PnL:";
audBad:.evidence.audit badRun;
show audBad`checks;
-1 "    overall pass = ",(string audBad`pass),"   (",(audBad`reason),")";
-1 "";

/ (3) the fail-safe precondition wrapper: reject the contaminated run; the gates never run.
.gov.register `hypoId`thesis`edgeSource`instruments`claimedRegimes`status!(
    `H_EVID;"evidence-audit demo";`riskPremium;enlist comm;enlist `backwardation;`proposed);
runnerCalls:0;
runner:{[from;to] runnerCalls+:1; sub:dts where dts within (from;to); ([] date:sub; pnl:(count sub)#0.0)};
-1 "(3) .evidence.gatedRun - the HARD precondition (mirrors carded gating):";
runnerCalls:0;
vBad:.evidence.gatedRun[`H_EVID;badRun;runner;`curveState];
-1 "    contaminated run -> verdict=",(string (first vBad)`verdict),", tradeable=",(string (first vBad)`tradeable),", gates run? ",string 0<runnerCalls;
runnerCalls:0;
@[{.evidence.gatedRun[`H_EVID;x;runner;`curveState]};cleanRun;{x}];
-1 "    clean run        -> audit passed, delegated to the gates (runner invoked? ",(string 0<runnerCalls),")";
-1 "";
-1 "The door (R9) prevents look-ahead, the replay (R12) records it, the audit (R13) PROVES it: a run";
-1 "that cannot be proven clean is REJECTED before the gates ever deflate a fabricated Sharpe.";

exit 0;
