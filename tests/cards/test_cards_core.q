\l core/init.q
/ ============================================================================
/ test_cards_core.q - model-card round-trip (.cards.get / .list / .forCapability).
/ Synthetic in-memory card set (no HDB). Research OS R5.
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};
caught:{[f;arg] @[f;arg;{[e] 1b}]~1b};

/ a synthetic card set (override the global the layer reads).
modelCards:(.cards.__empty[]) upsert `cardId`capabilityKind`capabilityName`version`intendedUse`assumptions`edgeSource`regimeApplicability`riskMemoryKey`govHypoId`owner`asOf!(
    `card_synSig;`signal;`synSig;`v1;"a synthetic signal";"some assumptions";`riskPremium;"trending regimes";`na;`;"tester";2026.05.31);
modelCards:modelCards upsert `cardId`capabilityKind`capabilityName`version`intendedUse`assumptions`edgeSource`regimeApplicability`riskMemoryKey`govHypoId`owner`asOf!(
    `card_synPricer;`model;`synPricer;`v1;"a synthetic pricer";"lognormal";`na;"any";`na;`;"tester";2026.05.31);

/ --- get resolves a card ---
c:.cards.get[`synSig];
chk[(c`capabilityName)~`synSig; "get must return the card for the capability"];
chk[(c`edgeSource)~`riskPremium; "the card's edgeSource must round-trip"];

/ --- list includes both ---
chk[2=count .cards.list[]; "list must include both synthetic cards"];
chk[`synSig in exec capabilityName from .cards.list[]; "list must include synSig"];

/ --- forCapability filters by kind + name ---
chk[1=count .cards.forCapability[`signal;`synSig]; "forCapability must resolve the (kind,name) card"];
chk[0=count .cards.forCapability[`model;`synSig]; "forCapability must not cross kinds"];

/ --- unknown capability -> signal ---
chk[caught[.cards.get;`nope]; "get of an unknown capability must signal"];

-1 "test_cards_core: PASS";
