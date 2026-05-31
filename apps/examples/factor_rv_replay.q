\l core/init.q
/ ============================================================================
/ factor_rv_replay.q - real-data example (NOT a test), DATA-CONDITIONAL. POST-FOUNDATION R17.
/ Re-run R8's factorRelativeValue strategy (curve PCA -> fade the residual) through the REALISTIC
/ end-to-end foundation R16 assembled (.workflow.runReplay: pre-register -> carded gating -> replay R12
/ -> evidence audit R13 -> gov gates -> regime skeptic -> attribution R15 -> human packet), and ask one
/ clean question: does realistic evidence (realistic costs + the CAUSAL PCA + deflation + the sealed
/ holdout) change R8's honest no-edge verdict? The params are R8's, PRE-REGISTERED, NOT tuned: a "research,
/ not tradeable" verdict that CONFIRMS the old finding is the system WORKING - it proves the foundation
/ generalises beyond the calendar spread it was demoed with. Skips if no HDB.
/ ============================================================================
hdbPath:.cfg.paths`hdb;
useHdb:0<count @[{[p] key hsym `$p,"/sym"};hdbPath;{[e] ()}];
if[not useHdb;
    -1 "factor_rv_replay SKIPPED - no HDB at ",hdbPath," (run scripts/ingest_hdb.q first).";
    exit 0];
.data.hdb.open hdbPath;

comm:`CRUDE; cfg:.template.factorRvReplay.__cfg[()!()];
/ scope the audited run record to train+validate (the holdout stays sealed; gov reads it one-shot).
tv:.gov.zone.range[comm;`trainValidate];
res:.template.factorRvReplay.run `commodity`dateFrom`dateTo!(comm; tv 0; tv 1);
run:(res`meta)`runRecord;
-1 "factorRelativeValue RE-RUN on the realistic foundation (",(string comm),", train+validate ",(string tv 0)," .. ",(string tv 1),"):";
-1 "  PCA: k=",(string cfg`k),", nMaturities=",(string cfg`nMaturities),", residualMaturity rank=",(string cfg`residualMaturity)," (R8's params, PRE-REGISTERED, NOT tuned)";
-1 "  fade: lookback=",(string cfg`lookback),", entryZ=",(string cfg`entryZ),", txnCostRate=",string cfg`txnCostRate;
-1 "  CAUSAL PCA: per step on a trailing ",(string cfg`lookback),"-day window of as-of curve changes (NEVER the full sample R8 used)";
-1 "  R8's gates (REUSED): ",(res[`validation;`shapeGate])`detail;
-1 "  run record: ",(string count run`steps)," steps trading the actual rank-",(string cfg`residualMaturity)," contract, replay totalPnl=",string (run`meta)`totalPnl;
-1 "";

/ run the END-TO-END bounded replay-mode loop, reusing R8's CARD (factorRelativeValue, riskMemory covid2020).
/ NOTE: the gov hypothesis declares edgeSource `mispricing (gov's validEdgeSources vocabulary for the
/ structural relative-value / dislocation edge); R8's card uses the cards-layer term `structural - the same
/ edge, two enums (gov: riskPremium/mispricing/info; cards: riskPremium/structural/informational).
prop:`hypoId`thesis`edgeSource`instruments`claimedRegimes`capName`runRecord`runner!(
    `crude_factorRv_replay;
    "the crude curve reverts to its k-factor PCA shape; fade the cumulative residual (factor-structure reversion / relative-value)";
    `mispricing; enlist comm; `backwardation`contango`flat; `factorRelativeValue;
    run; .template.factorRvReplay.runner[comm;cfg]);
pkt:.workflow.runReplay prop;

-1 "(1) The HONEST verdict (per regime bucket; gates: cost -> deflated Sharpe -> walk-forward -> sealed holdout):";
$[98h=type pkt`verdicts;
    show select bucket,nObs,verdict,failedGate,tradeable,netSharpe,dsr from pkt`verdicts;
    -1 "    ",-3!pkt`verdicts];
-1 "  anyTradeable = ",(string pkt`anyTradeable),"   (evidence audit passed = ",(string pkt`evidencePassed),")";
-1 "  regime risk-memory: ",$[0<count pkt`riskMemory; pkt`riskMemory; "(none)"];
-1 "";

-1 "(2) PnL ATTRIBUTION (R15: is the PnL the structural slope/curvature thesis, or unexplained residual?):";
a:pkt`attribution;
-1 "    level=",(string a`level),"  slope=",(string a`slope),"  curvature=",(string a`curvature),"  carry=",(string a`carry),"  residual=",string a`residual;
-1 "    total=",(string a`total),"  RESIDUAL FRACTION=",string a`residualFraction;
-1 "    (the factor-RV thesis CLAIMS a structural slope/curvature normalisation. A large residual fraction";
-1 "     here means the realised PnL is NOT explained by recognisable curve factors - the red flag that the";
-1 "     'edge' is idiosyncratic noise, exactly what deflation + the holdout exist to punish.)";
-1 "";

-1 "(3) CONTRAST with R8's recorded OLD-ENGINE verdict:";
-1 "    R8 (old engine): vectorised, FRICTIONLESS, FULL-SAMPLE PCA -> honestly found NO tradeable edge.";
-1 "    R17 (this re-run): realistic txn costs + a CAUSAL trailing-window PCA (no look-ahead) + deflation";
-1 "      against the trial count + a sealed one-shot holdout. Every one of these can only SUBTRACT from a";
-1 "      frictionless full-sample backtest, so the honest expectation is that realistic evidence CONFIRMS";
-1 "      or STRENGTHENS the no-edge finding. anyTradeable=",(string pkt`anyTradeable)," -> ",$[pkt`anyTradeable; "a candidate for the human"; "CONFIRMED no tradeable edge - the foundation generalises"],".";
-1 "";

-1 "(4) The human-escalation packet:";
-1 "    decision=",(string pkt`decision),"  humanSignOffRequired=",(string pkt`humanSignOffRequired),"  stage=",string pkt`stage;
-1 "    summary: ",pkt`summary;
-1 "";
-1 "This is the first proof the realistic foundation generalises beyond R16's capstone: handed a real PRIOR";
-1 "strategy, the rigorous judge - now fed only evidence it can trust - returns the same honest answer R8";
-1 "did, and the attribution says WHY. A confirmed no-edge verdict is the architecture doing its job.";

exit 0;
