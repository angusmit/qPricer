\l core/init.q
/ ============================================================================
/ test_gov_verdict_failsafe.q - the FAIL-SAFE verdict contract (Research OS R3b).
/ A gate FAILURE can only produce a NON-tradeable verdict (`reject` / `research`);
/ a tradeable verdict (`pass` / `regimeConditional`) requires passing ALL gates and
/ carries tradeable=1b. Downstream must read `tradeable`, never the label string.
/ Synthetic returns, no HDB.
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};

mkHypo:{[claimed] `hypoId`thesis`edgeSource`instruments`claimedRegimes`status!(
    `H;"t";`riskPremium;`SYN;claimed;`research)};
ev:{[hypo;rets;bucket;nT;V;cfg]
    .gov.evaluate `hypo`rets`bucket`axis`nTrials`varSR`cfg!(hypo;rets;bucket;`curveState;nT;V;cfg)};
cfg:.cfg.gov;

posBig:0.0008+0.0004*sin 0.3*til 500;   / strong, sign-stable edge
zeroR:-0.00002+0.0004*sin 0.3*til 500;  / ~zero drift -> fails cost
small:0.0004+0.002*sin 0.4*til 120;     / small sample -> fails deflation under N=20

/ collect a FAILED case from each failing gate + the two ALL-PASS cases.
failCost:ev[mkHypo `backwardation`contango;zeroR;`contango;1;0f;cfg];
failDefl:ev[mkHypo enlist `backwardation;small;`contango;20;0.01;cfg];
failThesis:ev[`hypoId`thesis`instruments`claimedRegimes`status!(`H;"x";`SYN;enlist `backwardation;`research);posBig;`contango;1;0f;cfg];
passClaimed:ev[mkHypo `backwardation`contango;posBig;`contango;1;0f;cfg];
passPostHoc:ev[mkHypo enlist `backwardation;posBig;`contango;1;0f;cfg];

failed:(failCost;failDefl;failThesis);
allPass:(passClaimed;passPostHoc);

/ --- EVERY failed-gate case: tradeable=0b AND label in {reject,research} ---
chk[not any {x`tradeable} each failed; "a FAILED gate must never be tradeable"];
chk[all {(x`verdict) in `reject`research} each failed; "a FAILED gate's verdict must be reject or research"];
chk[all {not `none~x`failedGate} each failed; "a failed case must name a failedGate"];

/ --- EVERY all-gates-pass case: tradeable=1b AND label in {pass,regimeConditional} ---
chk[all {x`tradeable} each allPass; "an all-gates-pass case must be tradeable=1b"];
chk[all {(x`verdict) in `pass`regimeConditional} each allPass; "a tradeable verdict must be pass or regimeConditional"];
chk[all {`none~x`failedGate} each allPass; "an all-gates-pass case has no failed gate"];

/ --- the specific fail-safe FIX: a post-hoc DEFLATION FAILURE is research (not the
/ tradeable regimeConditional), while a post-hoc ALL-PASS is the tradeable regimeConditional ---
chk[(failDefl`verdict)~`research; "post-hoc deflation FAILURE -> research (not regimeConditional)"];
chk[(passPostHoc`verdict)~`regimeConditional; "post-hoc ALL-PASS -> regimeConditional (tradeable)"];
chk[(failDefl`tradeable)<passPostHoc`tradeable; "the failure is non-tradeable; the all-pass is tradeable"];

-1 "test_gov_verdict_failsafe: PASS";
