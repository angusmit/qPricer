\l core/init.q
/ ============================================================================
/ test_cards_validation_status.q - .cards.validationStatus is DERIVED from the gov
/ ledger, never asserted. Synthetic card + synthetic gov ledger (no HDB). Research OS R5.
/ A logged, deflation-failing, no-holdout hypothesis -> NOT tradeable + deflationPass=0b;
/ a card with no govHypoId (or no ledger entry) -> `ungated`.
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};

/ --- synthetic cards: one linked to a gov hypothesis, one with no governance link ---
modelCards:(.cards.__empty[]) upsert `cardId`capabilityKind`capabilityName`version`intendedUse`assumptions`edgeSource`regimeApplicability`riskMemoryKey`govHypoId`owner`asOf!(
    `card_synMom;`signal;`synMom;`v1;"synthetic momentum";"trend persists";`riskPremium;"trends";`na;`synHypo;"tester";2026.05.31);
modelCards:modelCards upsert `cardId`capabilityKind`capabilityName`version`intendedUse`assumptions`edgeSource`regimeApplicability`riskMemoryKey`govHypoId`owner`asOf!(
    `card_synPricer;`model;`synPricer;`v1;"synthetic pricer";"lognormal";`na;"any";`na;`;"tester";2026.05.31);

/ --- reset the gov ledger and log a deflation-FAILING family for synHypo ---
.gov.hypoTbl:.gov.__emptyHypotheses[];
.gov.trialTbl:.gov.__emptyTrials[];
.gov.register `hypoId`thesis`edgeSource`instruments`claimedRegimes`status!(
    `synHypo;"synthetic";`riskPremium;`SYN;enlist `backwardation;`research);
/ three small-sample, modest-Sharpe trials -> recomputed DSR well below 0.95 (N=3, short n).
.gov.logTrial `hypoId`regimeBucket`nObs`netSharpe`skew`kurtosis!(`synHypo;`backwardation;120;0.6;0f;3f);
.gov.logTrial `hypoId`regimeBucket`nObs`netSharpe`skew`kurtosis!(`synHypo;`contango;90;0.9;0f;3f);
.gov.logTrial `hypoId`regimeBucket`nObs`netSharpe`skew`kurtosis!(`synHypo;`flat;70;1.4;0f;3f);

vstat:.cards.validationStatus[`synMom];
chk[(vstat`hypoId)~`synHypo; "validationStatus must link the card's govHypoId"];
chk[3=vstat`nTrials; "validationStatus must read the ledger trial count (N=3)"];
chk[(vstat`status)~`inResearch; "logged but no holdout look -> inResearch (NOT validated)"];
chk[not vstat`tradeable; "a deflation-failing, un-promoted capability must NOT be tradeable"];
chk[not vstat`deflationPass; "the recomputed headline DSR must be below the threshold (failed deflation)"];
chk[(vstat`headlineDsr)<0.95; "the derived DSR must reflect the small-sample penalty"];

/ --- a card with no governance link -> ungated (honest, never 'validated') ---
vp:.cards.validationStatus[`synPricer];
chk[(vp`status)~`ungated; "no govHypoId -> ungated"];
chk[not vp`tradeable; "an ungated capability is not tradeable"];

/ --- a HOLDOUT-PASSED hypothesis -> validated + tradeable (the only way to get there) ---
.gov.hypoTbl:update holdoutUsedAt:2026.05.31D00:00:00.0, holdoutVerdict:`pass, holdoutNetSharpe:1.2, holdoutDsr:0.97 from .gov.hypoTbl where hypoId=`synHypo;
vfin:.cards.validationStatus[`synMom];
chk[(vfin`status)~`validated; "a passed holdout look -> validated"];
chk[vfin`tradeable; "a passed holdout look -> tradeable"];
chk[vfin`holdoutUsed; "holdoutUsed must reflect the recorded one-shot look"];

-1 "test_cards_validation_status: PASS";
