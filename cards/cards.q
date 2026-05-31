/ cards/cards.q - model cards: the KNOWLEDGE plug-in (.cards.*) (v0.69, Research OS R5)
/ ----------------------------------------------------------------------------
/ ARCHITECTURE.md Part II (knowledge plug-ins) / 13-R5. A model card is a structured,
/ queryable record per capability that SYNTHESISES the rest of the system: what the
/ capability IS (the R2 contract it satisfies), what it assumes, which of the three edge
/ sources it claims (for signals/strategies), which regimes it's valid in, its linked
/ risk memory (R4), and - the part that makes a card honest rather than marketing - its
/ VALIDATION STATUS DERIVED FROM THE GOV LEDGER (R3/R3b), never asserted.
/ -
/ HIGH layer: a card sits ABOVE the things it documents. It READS (all legal downward):
/ the R2 registries (.model/.signal/.calibrator/.execution.fillModel.list + the strategy
/ registry), the gov ledger + hypotheses (.gov.trials / .gov.hypoTbl, recomputing the
/ deflated Sharpe via .gov.deflatedSharpe), and the R4 regime library (.regime.library).
/ It MUST NOT be imported by gov//regime//backtest/. Loads after gov/, before apps/.
/ Never opens the HDB at import; .cards.open uses `get`, not `\l`.
/ -
/ The honesty through-line continues: .cards.audit[] CAN FAIL (an undocumented registered
/ capability, a card missing a required section, a strategy/signal card with no named edge
/ source, an orphan card), and .cards.validationStatus CANNOT say "validated" when the
/ ledger says otherwise - it returns `ungated` when nothing is logged.
/ ----------------------------------------------------------------------------

/ Typed empty card table (the schema). The curated content lives in .cfg.cards; the
/ on-disk splay is `modelCards`. govHypoId links the card to its governance record (` if none).
.cards.__empty:{[]
    ([] cardId:`symbol$(); capabilityKind:`symbol$(); capabilityName:`symbol$(); version:`symbol$();
        intendedUse:(); assumptions:(); edgeSource:`symbol$(); regimeApplicability:();
        riskMemoryKey:`symbol$(); govHypoId:`symbol$(); owner:(); asOf:`date$())
 };

/ The active card set: the persisted `modelCards table if loaded, else the curated
/ .cfg.cards spec (so the knowledge plug-in works out-of-the-box; .cards.open / a build
/ override it). Tests set `modelCards directly.
.cards.cards:{[] $[`modelCards in key `.; modelCards; .cfg.cards]};

/ ── card lookup ──────────────────────────────────────────────────────────────

.cards.list:{[] select cardId,capabilityKind,capabilityName,version,edgeSource,asOf from .cards.cards[]};

/ NOTE: parameters are capName/capKind (NOT capabilityName/capabilityKind) - a qSQL
/ `where capabilityName=capabilityName` would compare the COLUMN to itself (always true).
.cards.get:{[capName]
    c:select from .cards.cards[] where capabilityName=capName;
    if[0=count c; '"cards.get: no card for ",string capName];
    first c
 };

.cards.forCapability:{[capKind;capName]
    select from .cards.cards[] where capabilityKind=capKind, capabilityName=capName
 };

/ ── validation status DERIVED from governance (never asserted) ───────────────

/ The capability's validation status read from the gov ledger via its card's govHypoId.
/ Recomputes the headline deflated Sharpe over the family's trials (same maths as
/ .gov.evaluate, reusing .gov.deflatedSharpe), and reads the holdout outcome from the
/ hypotheses registry. NOTHING logged (no govHypoId, or no trials) -> `ungated` (honest).
/ status: `ungated` (no record) / `inResearch` (logged, no holdout look yet) /
/ `validated` (holdout passed) / `rejected` (holdout failed). tradeable IFF holdout passed.
.cards.validationStatus:{[capName]
    base:`capabilityName`hypoId`status`nTrials`headlineNetSharpe`headlineDsr`deflationPass`holdoutUsed`tradeable!(
        capName;`;`ungated;0;0n;0n;0b;0b;0b);
    cs:select from .cards.cards[] where capabilityName=capName;
    hid:$[(0<count cs) and `govHypoId in cols cs; (first cs)`govHypoId; `];
    if[null hid; :base];
    base[`hypoId]:hid;                          / surface the link even when nothing is logged
    fam:.gov.trials hid;
    if[0=count fam; :base];
    annDays:.cfg.gov`annualizationDays;
    perPeriod:(fam`netSharpe)%sqrt annDays;
    nT:count fam;
    varSR:$[1<nT; var perPeriod; 0f];
    hi:fam first idesc fam`netSharpe;            / the headline (max net-Sharpe) trial
    srPP:(hi`netSharpe)%sqrt annDays;
    dsr:.gov.deflatedSharpe[srPP;hi`nObs;hi`skew;hi`kurtosis;nT;varSR];
    deflationPass:dsr>=.cfg.gov`dsrThreshold;
    hRows:select from .gov.hypoTbl where hypoId=hid;
    holdoutUsed:$[0<count hRows; not null (first hRows)`holdoutUsedAt; 0b];
    hv:$[0<count hRows; (first hRows)`holdoutVerdict; `];
    tradeable:holdoutUsed and hv=`pass;
    status:$[holdoutUsed; $[hv=`pass; `validated; `rejected]; `inResearch];
    `capabilityName`hypoId`status`nTrials`headlineNetSharpe`headlineDsr`deflationPass`holdoutUsed`tradeable!(
        capName;hid;status;nT;hi`netSharpe;dsr;deflationPass;holdoutUsed;tradeable)
 };

/ The linked regime-library risk memory for a card (reuses R4). Empty if the card has no
/ riskMemoryKey (`na` / `).
.cards.riskMemory:{[capName]
    cs:select from .cards.cards[] where capabilityName=capName;
    if[0=count cs; :.regime.library.__empty[]];
    rk:(first cs)`riskMemoryKey;
    if[(rk=`na) or null rk; :()];
    select episodeId,label,driversKey,riskMemory from .regime.library.episodes[] where driversKey=rk
 };

/ ── the audit (the honesty check that CAN FAIL) ──────────────────────────────

/ The registered capabilities across every kind. DYNAMIC (v0.71): it walks ALL kinds in
/ the R2 registry (`.registry.kinds[]`) - so a NEW kind (e.g. the R6 `template` kind, or
/ any future kind) can never be silently missed by the coverage audit - plus the separate
/ strategy registry. A new plug-in kind is visible to the audit the moment it is registered.
.cards.__registered:{[]
    rk:.registry.kinds[];
    (rk!{.registry.list x} each rk),(enlist `strategy)!enlist .strategy.registeredStrategies[]
 };

/ PURE audit: cross-check a registered-capabilities map (kind -> names) and a card set.
/ One row per registered capability (documented? required sections non-empty? edge named
/ for signal/strategy?) + one row per orphan card (capabilityName not registered).
.cards.__auditAgainst:{[registered;cards]
    validEdges:`riskPremium`structural`informational;
    rowFor:{[cards;validEdges;k;n]
        cs:select from cards where capabilityName=n;
        reasons:();
        if[0=count cs; reasons,:enlist "undocumented (no card)"];
        if[0<count cs;
            cc:first cs;
            if[0=count cc`intendedUse; reasons,:enlist "empty intendedUse"];
            if[0=count cc`assumptions; reasons,:enlist "empty assumptions"];
            if[0=count cc`regimeApplicability; reasons,:enlist "empty regimeApplicability"];
            if[(k in `signal`strategy) and not (cc`edgeSource) in validEdges;
                reasons,:enlist "no named edge source (need riskPremium/structural/informational)"]];
        `kind`name`pass`reason!(k;n;0=count reasons; $[0=count reasons; "ok"; "; " sv reasons])};
    regRows:raze {[rowFor;cards;validEdges;k;ns] rowFor[cards;validEdges;k;] each ns}[rowFor;cards;validEdges]'[key registered; value registered];
    allNames:raze value registered;
    orphans:select from cards where not capabilityName in allNames;
    orphanRows:$[0=count orphans; ();
        {[oc] `kind`name`pass`reason!(oc`capabilityKind;oc`capabilityName;0b;"orphan card - capabilityName not registered")} each orphans];
    regRows,orphanRows
 };

/ Audit the live cards against the live registries.
.cards.audit:{[] .cards.__auditAgainst[.cards.__registered[]; .cards.cards[]]};

/ ── connective wiring (v0.71): carded gating ────────────────────────────────

/ Is a capability ready to be gated? -> `ready`reason. Ready iff a model card EXISTS for it
/ AND its failure-mode / risk-memory field is populated (riskMemoryKey is a real link, not
/ `na / empty - it resolves to the regime library's documented failure modes). Honest gate:
/ you cannot gate a capability you have not documented the failure modes of.
.cards.gateReady:{[capName]
    cs:select from .cards.cards[] where capabilityName=capName;
    if[0=count cs; :`ready`reason!(0b;"no model card for ",string capName)];
    rk:(first cs)`riskMemoryKey;
    populated:(not null rk) and not rk=`na;
    $[populated;
        `ready`reason!(1b;"carded with a populated failure-mode (risk-memory) link: ",string rk);
        `ready`reason!(0b;"card has no populated failure-mode field (riskMemoryKey is `na / empty)")]
 };

/ Carded gating ENTRY: realises "Gate 0 requires a card with a populated failure-mode field"
/ at the BOUNDARY (cards -> gov is a legal DOWNWARD call; gov never imports cards). Prechecks
/ .cards.gateReady; if NOT ready it REFUSES - returns an `undocumented verdict (tradeable=0b)
/ and does NOT run the gates / does NOT touch the holdout. If ready it delegates to
/ .gov.runFull and returns its verdict (now carrying the regime risk-memory annotation).
.cards.gatedRun:{[hypoId;capName;runner;axis]
    gr:.cards.gateReady capName;
    if[not gr`ready;
        :enlist `bucket`nObs`verdict`failedGate`tradeable`netSharpe`dsr`postHoc`holdoutPassed`reason`riskMemory!(
            `all;0;`undocumented;`card;0b;0n;0n;0b;0b;("REFUSED before gating: ",gr`reason);"")];
    .gov.runFull[hypoId;runner;axis]
 };

/ ── HDB persistence (touches real gitignored data; build script / demo only) ─

/ Write `modelCards` from the curated .cfg.cards spec (mirrors scripts/build_regime_library.q:
/ .Q.en + set, via get/`set`, no \l). Idempotent.
.cards.buildTable:{[hdbPath]
    t:.Q.en[hsym `$hdbPath; .cfg.cards];
    (hsym `$hdbPath,"/modelCards/") set t;
    count t
 };

/ Load the persisted modelCards (+ sym domain) with `get` - NOT \l (no cwd change).
.cards.open:{[hdbPath]
    symPath:hsym `$hdbPath,"/sym";
    if[0<count key symPath; sym::get symPath];
    mp:hsym `$hdbPath,"/modelCards/";
    if[0=count key mp; '"cards.open: no modelCards at ",hdbPath," (run scripts/build_model_cards.q)"];
    modelCards::get mp;
    hdbPath
 };

-1 "cards.q loaded - .cards.* model-card knowledge plug-in ready (",(string count .cards.cards[])," curated cards)";
