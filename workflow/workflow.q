/ workflow/workflow.q - the single-process research loop (.workflow.*) (v0.72, Research OS R7)
/ ----------------------------------------------------------------------------
/ ARCHITECTURE.md Part II 11.7 / 13-R7. A THIN HIGH layer that COMPOSES the already-built
/ functions into the bounded research loop the six agent roles (agents/*.md) enact:
/   researcher proposes  -> curator ensures carded  -> validator gates (one-shot holdout)
/   -> skeptic annotates  -> logger records  -> lead packages + ESCALATES to the human.
/ It REUSES .gov.* / .cards.* / .template.* / .regime.* - it does NOT reimplement anything,
/ adds NO compute path, and opens NO HDB at import. It composes DOWNWARD only (it is the top
/ functional layer; loads after cards/, before apps/).
/ -
/ BOUNDED BY CONSTRUCTION (the defining R7 principle): the loop has NO path that allocates
/ capital, trades, or marks something tradeable without ALL gates passing. qFDM is a BATCH
/ research system - there is no order execution. The workflow PACKAGES the honest case for a
/ HUMAN; the human holds every go/no-go. The escalation packet has no deploy/allocate field.
/ -
/ The carded-gating + skeptic-annotation + trial-logging are ALREADY composed inside
/ .cards.gatedRun -> .gov.runFull (the v0.71 wiring); .workflow.run sequences the roles
/ around that and assembles the escalation packet.
/ ----------------------------------------------------------------------------

/ Assemble the human-escalation packet. carded/gatesRun reflect whether the capability was
/ card-ready and the gates were run; verdicts is the per-bucket gate result (a table) or an
/ empty dict on refusal. anyTradeable is DERIVED (true iff some bucket passed ALL gates -
/ the R3b fail-safe) - it is a FACT, not a decision. The packet ALWAYS defers to the human
/ (decision=`escalateToHuman, humanSignOffRequired=1b) and starts staging at `paper; it
/ carries NO deploy / allocate / autoTrade field by construction.
.workflow.__packet:{[hypoId;capName;carded;verdicts;reason]
    isTbl:(98h=type verdicts) and 0<count verdicts;
    anyTradeable:$[isTbl; any verdicts`tradeable; 0b];
    rmAll:$[isTbl; verdicts`riskMemory; ()];
    nonEmpty:rmAll where 0<count each rmAll;
    riskMemory:$[0<count nonEmpty; first nonEmpty; ""];
    summary:$[not carded;
        "REFUSED at the carded step: ",reason," - document the capability's failure modes before it can be gated";
        $[anyTradeable;
            "cleared ALL gates in >=1 regime bucket -> a CANDIDATE for the human (NOT auto-deployed); paper stage first";
            "did NOT clear all gates -> research only; no tradeable claim"]];
    `hypoId`capability`carded`gatesRun`anyTradeable`verdicts`riskMemory`decision`humanSignOffRequired`stage`summary!(
        hypoId;capName;carded;carded;anyTradeable;verdicts;riskMemory;`escalateToHuman;1b;`paper;summary)
 };

/ Run the bounded research loop for a proposal and return the human-escalation packet.
/ proposal (a dict): hypoId, thesis, edgeSource (riskPremium/structural/informational/...),
/ instruments (the commodity/instruments), claimedRegimes, capName (the capability to gate -
/ matching a model card + a registry name), runner (a runner[from;to]->(date;pnl) callback),
/ optional axis (the regime axis, default `curveState).
.workflow.run:{[proposal]
    hypoId:proposal`hypoId;
    capName:proposal`capName;
    axis:$[`axis in key proposal; proposal`axis; `curveState];
    / 1. RESEARCHER: pre-register the hypothesis (thesis + edge source + claimed regimes),
    / BEFORE looking at any outcome. Idempotent.
    .gov.register `hypoId`thesis`edgeSource`instruments`claimedRegimes`status!(
        hypoId;
        proposal`thesis;
        proposal`edgeSource;
        proposal`instruments;
        proposal`claimedRegimes;
        `research);
    / 2. CURATOR: the capability must be carded WITH a populated failure-mode field before it
    / can be gated. If not card-ready, REFUSE here - the gates are NOT run, the holdout is
    / NEVER touched.
    gr:.cards.gateReady capName;
    if[not gr`ready;
        :.workflow.__packet[hypoId;capName;0b;()!();gr`reason]];
    / 3+4+5. VALIDATOR gate (one-shot holdout) + SKEPTIC riskMemory annotation + LOGGER trial
    / record - all already composed by .cards.gatedRun -> .gov.runFull.
    verdicts:.cards.gatedRun[hypoId;capName;proposal`runner;axis];
    / 6. LEAD: package the honest case and escalate to the human.
    .workflow.__packet[hypoId;capName;1b;verdicts;"gated"]
 };

-1 "workflow.q loaded - .workflow.* bounded research loop ready (composes gov/cards/template/regime)";
