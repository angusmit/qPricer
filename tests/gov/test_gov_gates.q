\l core/init.q
/ ============================================================================
/ test_gov_gates.q - the ordered GATE CASCADE (.gov.evaluate). Synthetic
/ (hypothesis, net-return vector, trial-count + Sharpe-variance) -> KNOWN
/ verdicts. Asserts stop-at-first-failure, the verdict-record fields, and the
/ Gate-0 post-hoc flag (tested bucket not in claimedRegimes). Research OS R3.
/ Deterministic sine-based returns -> reproducible (no RNG). No HDB.
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};

mkHypo:{[claimed] `hypoId`thesis`edgeSource`instruments`claimedRegimes`status!(
    `H;"crude TS momentum earns a positive excess return";`riskPremium;`CRUDE;claimed;`research)};
ev:{[hypo;rets;bucket;nTrials;varSR;cfg]
    .gov.evaluate `hypo`rets`bucket`axis`nTrials`varSR`cfg!(hypo;rets;bucket;`curveState;nTrials;varSR;cfg)};

cfg:.cfg.gov;
/ a Gate-3 override: an unreachable OOS Sharpe floor isolates the walk-forward gate.
wfCfg:@[cfg;`walkForward;:;(cfg`walkForward),(enlist `oosSharpeFloor)!enlist 100f];

/ deterministic return vectors (always-positive drift / ~zero drift / modest short sample).
posBig:0.0008+0.0004*sin 0.3*til 500;   / strong, sign-stable positive edge, n=500
zeroR:-0.00002+0.0004*sin 0.3*til 500;  / ~zero / slightly negative drift, n=500
small:0.0004+0.002*sin 0.4*til 120;     / positive but small sample, n=120

hClaimContango:mkHypo `backwardation`contango;   / contango IS claimed -> not post-hoc
hClaimBackOnly:mkHypo enlist `backwardation;      / contango NOT claimed -> post-hoc
hNoEdge:`hypoId`thesis`instruments`claimedRegimes`status!(`H;"x";`CRUDE;enlist `backwardation;`research);

/ --- 1. passes ALL gates ---
v1:ev[hClaimContango;posBig;`contango;1;0f;cfg];
chk[(v1`verdict)~`pass; "case1 verdict must be pass"];
chk[(v1`failedGate)~`none; "case1 no failed gate"];
chk[4=count v1`passedGates; "case1 must pass all four gates"];
chk[not v1`postHoc; "case1 contango is claimed -> not post-hoc"];

/ --- 2. fails Gate 1 (cost) -> stop at first failure (only Gate 0 passed) ---
v2:ev[hClaimContango;zeroR;`contango;1;0f;cfg];
chk[(v2`verdict)~`reject; "case2 verdict must be reject"];
chk[(v2`failedGate)~`cost; "case2 must fail the cost gate"];
chk[1=count v2`passedGates; "case2 stop-at-first-failure: only thesis gate passed"];
chk[(v2`passedGates)~enlist `thesis; "case2 passedGates must be just (thesis)"];

/ --- 3. fails Gate 2 (deflated Sharpe), bucket CLAIMED -> research ---
v3:ev[hClaimContango;small;`contango;20;0.01;cfg];
chk[(v3`verdict)~`research; "case3 (claimed) deflation fail -> research"];
chk[(v3`failedGate)~`deflatedSharpe; "case3 must fail the deflated-Sharpe gate"];
chk[2=count v3`passedGates; "case3 thesis+cost passed before deflation"];
chk[(v3`dsr)<0.95; "case3 DSR must be below the 0.95 threshold"];
chk[not v3`postHoc; "case3 contango is claimed -> not post-hoc"];

/ --- 4. fails Gate 2, bucket NOT claimed -> regimeConditional + post-hoc flag ---
v4:ev[hClaimBackOnly;small;`contango;20;0.01;cfg];
chk[(v4`verdict)~`regimeConditional; "case4 (post-hoc) deflation fail -> regimeConditional"];
chk[v4`postHoc; "case4 contango not in claimedRegimes -> post-hoc flagged"];
chk[(v4`failedGate)~`deflatedSharpe; "case4 must fail the deflated-Sharpe gate"];

/ --- 5. fails Gate 3 (walk-forward) under an unreachable OOS floor ---
v5:ev[hClaimContango;posBig;`contango;1;0f;wfCfg];
chk[(v5`verdict)~`research; "case5 walk-forward fail -> research"];
chk[(v5`failedGate)~`walkForward; "case5 must fail the walk-forward gate"];
chk[3=count v5`passedGates; "case5 thesis+cost+deflatedSharpe passed before walk-forward"];

/ --- 6. fails Gate 0 (thesis): no valid edgeSource -> reject, nothing else evaluated ---
v6:ev[hNoEdge;posBig;`contango;1;0f;cfg];
chk[(v6`verdict)~`reject; "case6 missing edgeSource -> reject"];
chk[(v6`failedGate)~`thesis; "case6 must fail the thesis gate"];
chk[0=count v6`passedGates; "case6 stop-at-first-failure: no gate passed"];

-1 "test_gov_gates: PASS";
