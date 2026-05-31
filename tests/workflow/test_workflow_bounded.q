\l core/init.q
/ ============================================================================
/ test_workflow_bounded.q - the BOUNDEDNESS invariants (Research OS R7). The loop NEVER
/ returns anyTradeable=1b unless a bucket passed ALL gates (the R3b fail-safe holds end-to-
/ end), and the escalation packet ALWAYS defers the decision to the human - there is no
/ auto-allocation / auto-deploy field, ever. Synthetic; .gov.runFull stubbed (no HDB).
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};

modelCards:([] cardId:enlist `c_r; capabilityKind:enlist `signal; capabilityName:enlist `readyCap;
    version:enlist `v1; intendedUse:enlist "u"; assumptions:enlist "a"; edgeSource:enlist `riskPremium;
    regimeApplicability:enlist "t"; riskMemoryKey:enlist `energyShock2022; govHypoId:enlist `na;
    owner:enlist "t"; asOf:enlist 2026.05.31);
prop:`hypoId`thesis`edgeSource`instruments`claimedRegimes`capName`runner`axis!(
    `hB;"t";`riskPremium;`CRUDE;enlist `backwardation;`readyCap;
    {[from;to] ([] date:enlist 2020.01.01; pnl:enlist 0f)};`curveState);

origRunFull:.gov.runFull;

/ --- case 1: NO bucket passes all gates -> anyTradeable MUST be 0b ---
.gov.hypoTbl:.gov.__emptyHypotheses[]; .gov.trialTbl:.gov.__emptyTrials[];
.gov.runFull:{[h;r;a] ([] bucket:`b1`b2; nObs:100 50; verdict:`reject`research; failedGate:`cost`deflatedSharpe;
    tradeable:00b; netSharpe:-0.2 0.3; dsr:0n,0.4; postHoc:00b; holdoutPassed:00b; reason:("x";"y"); riskMemory:("";""))};
p0:.workflow.run prop;
chk[not p0`anyTradeable; "no gate-passing bucket -> anyTradeable MUST be 0b"];
chk[(p0`decision)~`escalateToHuman; "the packet defers to the human"];
chk[p0`humanSignOffRequired; "human sign-off is always required"];

/ --- case 2: a bucket passes all gates -> anyTradeable=1b, but STILL defers to the human ---
.gov.hypoTbl:.gov.__emptyHypotheses[]; .gov.trialTbl:.gov.__emptyTrials[];
.gov.runFull:{[h;r;a] ([] bucket:enlist `b1; nObs:enlist 300; verdict:enlist `pass; failedGate:enlist `none;
    tradeable:enlist 1b; netSharpe:enlist 1.8; dsr:enlist 0.98; postHoc:enlist 0b; holdoutPassed:enlist 1b;
    reason:enlist "ok"; riskMemory:enlist "resembles 2022: late entries round-tripped")};
p1:.workflow.run prop;
chk[p1`anyTradeable; "a bucket passing ALL gates -> anyTradeable=1b"];
chk[(p1`decision)~`escalateToHuman; "even when tradeable, the decision STILL defers to the human"];
chk[p1`humanSignOffRequired; "tradeable still requires human sign-off"];
chk[(p1`stage)~`paper; "capital staging starts at paper (advanced only by the human)"];

/ --- the packet has NO autonomous deploy / allocate / trade field, in EITHER case ---
forbidden:`deploy`allocate`autoTrade`autoDeploy`capital`order`execute;
chk[not any forbidden in key p0; "refused/non-tradeable packet carries NO auto-deploy/allocate field"];
chk[not any forbidden in key p1; "tradeable packet carries NO auto-deploy/allocate field (the human deploys)"];

.gov.runFull:origRunFull;
-1 "test_workflow_bounded: PASS";
