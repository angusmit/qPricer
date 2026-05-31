\l core/init.q
/ ============================================================================
/ test_cards_audit.q - the audit that CAN FAIL (.cards.__auditAgainst, the pure core
/ of .cards.audit). Synthetic registered-capability map + synthetic card set (no HDB).
/ Research OS R5. The audit must BITE: it flags (a) an undocumented registered
/ capability, (b) a card missing a required section, (c) a strategy card with no named
/ edge source, and (d) an orphan card - while passing fully-documented capabilities.
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};

/ Build the synthetic card set as ONE table literal (row-by-row upsert of an EMPTY-string
/ section would hit a column-type mismatch). missSec has empty assumptions; naEdge is a
/ strategy with edgeSource `na; ghostCap is an orphan (not in the registered map below).
cards:([]
    cardId:`c_docStrat`c_missSec`c_naEdge`c_pricer`c_ghost;
    capabilityKind:`strategy`strategy`strategy`model`strategy;
    capabilityName:`docStrat`missSec`naEdge`pricer1`ghostCap;
    version:5#`v1;
    intendedUse:("use";"use";"use";"use";"use");
    assumptions:("assume";"";"assume";"assume";"assume");
    edgeSource:`riskPremium`riskPremium`na`na`riskPremium;
    regimeApplicability:("trends";"trends";"trends";"any";"trends");
    riskMemoryKey:5#`na;
    govHypoId:5#`na;
    owner:5#enlist "t";
    asOf:5#2026.05.31);

/ a synthetic registered map (do NOT depend on the live registries).
registered:`model`strategy!((enlist `pricer1); `docStrat`missSec`naEdge`undocStrat);

a:.cards.__auditAgainst[registered;cards];
passOf:{[a;nm] first exec pass from a where name=nm};

/ --- fully-documented capabilities PASS ---
chk[passOf[a;`docStrat]; "a fully-documented strategy card must pass"];
chk[passOf[a;`pricer1]; "a model card with edgeSource `na must pass (na is fine for a pricer)"];

/ --- the audit BITES on all four failure modes ---
chk[not passOf[a;`undocStrat]; "(a) an undocumented registered capability must be FLAGGED"];
chk[not passOf[a;`missSec]; "(b) a card missing a required section must be FLAGGED"];
chk[not passOf[a;`naEdge]; "(c) a strategy card with no named edge source must be FLAGGED"];
chk[not passOf[a;`ghostCap]; "(d) an orphan card (capabilityName not registered) must be FLAGGED"];

/ --- the failure reasons are specific (not a blanket fail) ---
chk[(first exec reason from a where name=`undocStrat) like "*no card*"; "undocumented reason must say so"];
chk[(first exec reason from a where name=`naEdge) like "*edge*"; "no-edge reason must mention the edge source"];

/ --- coverage: 5 registered (1 model + 4 strategy) + 1 orphan = 6 rows; docStrat + pricer1
/ pass; missSec/naEdge/undocStrat + the ghost orphan fail = 4 fails ---
chk[6=count a; "audit must cover every registered capability + every orphan card"];
chk[2=sum a`pass; "exactly the two well-formed capabilities pass"];
chk[4=sum not a`pass; "the audit must flag all four failure modes"];

-1 "test_cards_audit: PASS";
