\l core/init.q
/ ============================================================================
/ test_gov_ledger.q - hypothesis registry + APPEND-ONLY trials ledger
/ (.gov.register / .gov.logTrial / .gov.trials / .gov.nTrials). Synthetic
/ in-memory tables only (no CSV, no HDB read). Research OS R3 (gov layer).
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};

/ start from EMPTY in-memory tables (test independence within the shared suite process).
.gov.hypoTbl:.gov.__emptyHypotheses[];
.gov.trialTbl:.gov.__emptyTrials[];

/ --- register a hypothesis; the id comes back ---
hid:.gov.register `hypoId`thesis`edgeSource`instruments`claimedRegimes`status!(
    `mom_crude;"crude TS momentum earns a positive excess return";`riskPremium;`CRUDE;enlist `backwardation;`research);
chk[hid~`mom_crude; "register must return the hypoId"];
chk[1=count .gov.hypoTbl; "one hypothesis registered"];
chk[(.gov.hypo[`mom_crude]`edgeSource)~`riskPremium; "hypo lookup must round-trip the edgeSource"];

/ --- register is idempotent by hypoId (re-register replaces, does not duplicate) ---
.gov.register `hypoId`thesis`edgeSource`claimedRegimes`status!(
    `mom_crude;"crude TS momentum (revised)";`riskPremium;enlist `backwardation;`research);
chk[1=count .gov.hypoTbl; "re-register must not duplicate the hypothesis row"];

/ --- log a trial; the ledger grows by one ---
.gov.logTrial `hypoId`regimeBucket`nObs`netSharpe!(`mom_crude;`backwardation;1460;0.5);
chk[1=.gov.nTrials `mom_crude; "one trial logged"];

/ --- THE append-only contract: a SECOND log ADDS a row, never overwrites ---
.gov.logTrial `hypoId`regimeBucket`nObs`netSharpe!(`mom_crude;`contango;249;0.9);
chk[2=.gov.nTrials `mom_crude; "second logTrial must APPEND (count -> 2)"];
fam:.gov.trials `mom_crude;
chk[2=count fam; "trials[hypoId] must return both rows"];
chk[(asc fam`netSharpe)~asc 0.5 0.9; "both trials' netSharpe values preserved (no overwrite)"];
chk[(asc fam`regimeBucket)~asc `backwardation`contango; "both regime buckets preserved"];

/ --- a different family has zero trials ---
chk[0=.gov.nTrials `other_hypo; "an unlogged family has nTrials 0"];

/ --- a third log keeps appending (monotone N) ---
.gov.logTrial `hypoId`regimeBucket`nObs`netSharpe!(`mom_crude;`flat;157;2.5);
chk[3=.gov.nTrials `mom_crude; "third logTrial -> N=3 (the multiple-testing denominator grows)"];

-1 "test_gov_ledger: PASS";
