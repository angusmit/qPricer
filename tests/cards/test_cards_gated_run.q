\l core/init.q
/ ============================================================================
/ test_cards_gated_run.q - carded gating (.cards.gateReady / .cards.gatedRun). v0.71 wiring.
/ "Gate 0 requires a card with a populated failure-mode field", realised at the cards->gov
/ boundary (gov never imports cards). Synthetic cards; the delegate path stubs .gov.runFull
/ so no HDB is needed. Asserts the REFUSE path does NOT enter the gov gates.
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};

modelCards:([]
    cardId:`c_ready`c_bare;
    capabilityKind:`signal`signal;
    capabilityName:`readyCap`bareCap;
    version:`v1`v1;
    intendedUse:("u";"u");
    assumptions:("a";"a");
    edgeSource:`riskPremium`riskPremium;
    regimeApplicability:("trends";"trends");
    riskMemoryKey:`someEpisode`na;          / readyCap links failure modes; bareCap does not
    govHypoId:`na`na;
    owner:2#enlist "t";
    asOf:2#2026.05.31);

/ --- gateReady ---
chk[(.cards.gateReady `readyCap)`ready; "a card with a populated failure-mode link is gate-ready"];
chk[not (.cards.gateReady `bareCap)`ready; "a card with riskMemoryKey `na is NOT gate-ready"];
chk[not (.cards.gateReady `noCardCap)`ready; "an uncarded capability is NOT gate-ready"];

/ --- gatedRun REFUSES an undocumented capability WITHOUT entering the gov gates ---
.gov.hypoTbl:.gov.__emptyHypotheses[];
.gov.trialTbl:.gov.__emptyTrials[];
ref:.cards.gatedRun[`hX;`bareCap;{[from;to] ([] date:enlist 2020.01.01; pnl:enlist 0f)};`curveState];
chk[(first ref`verdict)~`undocumented; "gatedRun must REFUSE an undocumented capability (verdict undocumented)"];
chk[(first ref`failedGate)~`card; "the refusal names the card gate"];
chk[not first ref`tradeable; "a refused run is never tradeable"];
chk[0=count .gov.trialTbl; "REFUSED -> the gov gates were NOT entered (no trial logged)"];

/ --- gatedRun DELEGATES to .gov.runFull when ready (stub runFull so no HDB is needed) ---
origRunFull:.gov.runFull;
.gov.runFull:{[hid;runner;axis] enlist `bucket`verdict`tradeable`stubCalled!(`all;`pass;1b;1b)};
out:.cards.gatedRun[`hY;`readyCap;{[from;to] ([] date:enlist 2020.01.01; pnl:enlist 0f)};`curveState];
chk[(first out)`stubCalled; "gatedRun must DELEGATE to .gov.runFull when the capability is gate-ready"];
.gov.runFull:origRunFull;                  / restore the real runFull

-1 "test_cards_gated_run: PASS";
