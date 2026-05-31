\l core/init.q
/ ============================================================================
/ test_workflow_loop.q - the single-process research loop (.workflow.run). Research OS R7.
/ A CARDED proposal runs the full loop and returns a well-formed escalation packet; an
/ UNCARDED proposal is REFUSED at the carded step (gates NOT run, holdout untouched).
/ Synthetic cards; .gov.runFull is stubbed so no HDB is needed.
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};

modelCards:([] cardId:`c_r`c_b; capabilityKind:`signal`signal; capabilityName:`readyCap`bareCap;
    version:`v1`v1; intendedUse:("u";"u"); assumptions:("a";"a"); edgeSource:`riskPremium`riskPremium;
    regimeApplicability:("t";"t"); riskMemoryKey:`energyShock2022`na; govHypoId:`na`na; owner:2#enlist "t"; asOf:2#2026.05.31);

mkProp:{[cap] `hypoId`thesis`edgeSource`instruments`claimedRegimes`capName`runner`axis!(
    `hW;"synthetic thesis";`riskPremium;`CRUDE;enlist `backwardation;cap;
    {[from;to] ([] date:enlist 2020.01.01; pnl:enlist 0f)};`curveState)};

.gov.hypoTbl:.gov.__emptyHypotheses[]; .gov.trialTbl:.gov.__emptyTrials[];

/ --- UNCARDED proposal -> REFUSED at the carded step, gates NOT run ---
ref:.workflow.run mkProp `bareCap;
chk[not ref`carded; "an uncarded capability is not carded"];
chk[not ref`gatesRun; "REFUSED -> the gates were NOT run"];
chk[not ref`anyTradeable; "a refused proposal is never tradeable"];
chk[0=count .gov.trialTbl; "REFUSED -> no trial logged (the gov gates were not entered, holdout untouched)"];
chk[(ref`decision)~`escalateToHuman; "even a refusal escalates to the human"];

/ --- CARDED proposal -> full loop -> well-formed escalation packet (stub runFull) ---
origRunFull:.gov.runFull;
.gov.runFull:{[h;r;a] ([] bucket:`b1`b2; nObs:100 50; verdict:`pass`research; failedGate:`none`deflatedSharpe;
    tradeable:10b; netSharpe:1.5 0.2; dsr:0.97 0.4; postHoc:00b; holdoutPassed:10b;
    reason:("ok";"x"); riskMemory:("resembles 2020: storage blew up";""))};
pk:.workflow.run mkProp `readyCap;
.gov.runFull:origRunFull;

chk[pk`carded; "a carded capability passes the carded step"];
chk[pk`gatesRun; "the gates were run"];
chk[(`hypoId`capability`carded`gatesRun`anyTradeable`verdicts`riskMemory`decision`humanSignOffRequired`stage`summary)~key pk;
    "the escalation packet has the expected fields"];
chk[98h=type pk`verdicts; "the packet carries the per-bucket verdict table"];
chk[pk`anyTradeable; "a gate-passing bucket -> anyTradeable=1b"];
chk[(pk`riskMemory)~"resembles 2020: storage blew up"; "the skeptic's riskMemory annotation is surfaced"];
chk[(pk`hypoId)~`hW; "the packet records the hypothesis id"];

-1 "test_workflow_loop: PASS";
