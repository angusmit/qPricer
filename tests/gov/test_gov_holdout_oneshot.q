\l core/init.q
/ ============================================================================
/ test_gov_holdout_oneshot.q - the sealed holdout gate is ONE-SHOT (.gov.holdoutGate).
/ Synthetic in-memory `futures` (so .gov.zone.range / .gov.holdout.read resolve a
/ holdout date range with no HDB) + a counting runner. Research OS R3b (gov layer).
/ Asserts: the FIRST look runs the runner, records holdoutUsedAt + a ledger row with
/ dataZone=`holdout; the SECOND look returns the SAME recorded verdict WITHOUT
/ recomputing (the runner counter does not increment, the timestamp is unchanged);
/ and the runner is only ever called on the HOLDOUT range (never a trainValidate date).
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};

/ synthetic futures for commodity SYN: 100 contiguous dates (no HDB file needed).
futures:([] commodity:`SYN; date:2022.01.01+til 100; settle:100f; volume:1000);

.gov.hypoTbl:.gov.__emptyHypotheses[];
.gov.trialTbl:.gov.__emptyTrials[];
.gov.register `hypoId`thesis`edgeSource`instruments`claimedRegimes`status!(
    `HS;"synthetic strategy that clears gates 0-3";`riskPremium;`SYN;enlist `backwardation;`research);

/ the expected holdout range (so we can assert the runner only ever saw holdout dates).
hRange:.gov.holdout.read[`SYN];

/ a counting runner: records its call count + every (from;to) range it was asked for,
/ and returns a positive-Sharpe net PnL series over the requested range.
calls:0;
seenRanges:();
runner:{[from;to]
    calls+::1;
    seenRanges,::enlist (from;to);
    nD:1+`int$to-from;
    ([] date:from+til nD; pnl:0.001+0.0004*sin 0.3*til nD)};

/ --- FIRST look: fresh, records the seal ---
r1:.gov.holdoutGate[`HS;runner];
chk[r1`fresh; "first holdoutGate call must be a FRESH look"];
chk[r1`passed; "first look: the positive-Sharpe synthetic must pass the holdout bar"];
chk[1=calls; "first look must call the runner exactly once"];
ts1:(.gov.hypo `HS)`holdoutUsedAt;
chk[not null ts1; "first look must record holdoutUsedAt"];
chk[1=count select from .gov.trialTbl where dataZone=`holdout; "first look must append ONE holdout ledger row"];
chk[(hRange)~first seenRanges; "the runner must be called on the HOLDOUT range"];

/ --- SECOND look: the seal - recorded verdict, NO recomputation ---
r2:.gov.holdoutGate[`HS;runner];
chk[not r2`fresh; "second holdoutGate call must be RECORDED (not fresh)"];
chk[(r2`passed)~r1`passed; "second look must return the SAME recorded verdict"];
chk[1=calls; "second look must NOT call the runner again (one look per hypothesis, ever)"];
chk[ts1~(.gov.hypo `HS)`holdoutUsedAt; "second look must NOT change the recorded timestamp"];
chk[1=count select from .gov.trialTbl where dataZone=`holdout; "second look must NOT append another holdout row"];

/ --- the runner was NEVER asked for a trainValidate date (the seal protects 0-3) ---
tvTo:(.gov.zone.range[`SYN;`trainValidate]) 1;
chk[all (hRange 0) <= seenRanges[;0]; "the holdout runner must never start before the holdout range"];
chk[all tvTo < seenRanges[;0]; "the holdout runner must never touch a trainValidate date"];

-1 "test_gov_holdout_oneshot: PASS";
