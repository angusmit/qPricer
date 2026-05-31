\l core/init.q
/ ============================================================================
/ factor_relative_value.q - real-data example (NOT a test), DATA-CONDITIONAL. Research OS R8.
/ The FIRST new research shape on the completed spine: factor-structure relative value.
/ (1) run .factor.* on the real crude curve - print the level/slope/curvature loadings +
/ explained variance (a sanity check the factors look like a term structure); (2) instantiate
/ the factorRelativeValue template - show its factor-stability + residual-stationarity gates;
/ (3) if those pass, run the PnL through the FULL bounded workflow (carded gating -> the gov
/ cost/deflation/walk-forward/holdout cascade -> the regime risk-memory skeptic) for an honest
/ verdict + the human-escalation packet. Do NOT tune k or the lookback. Skips if no HDB. exit 0;.
/ ============================================================================
hdbPath:.cfg.paths`hdb;
useHdb:0<count @[{[p] key hsym `$p,"/sym"};hdbPath;{[e] ()}];
if[not useHdb;
    -1 "factor_relative_value SKIPPED - no HDB at ",hdbPath," (run scripts/ingest_hdb.q first).";
    exit 0];
.data.hdb.open hdbPath;
if[0<count key hsym `$hdbPath,"/regimeEpisodes/"; .regime.library.open hdbPath];

dates:.data.hdb.dates `CRUDE;
clHist:.data.hdb.curveHistory[`CRUDE;dates];

/ (1) the factor decomposition of the real crude curve.
d:.factor.decompose[`CRUDE;dates;.cfg.factor];
-1 "Curve PCA (CRUDE, ",(string .cfg.factor`nMaturities)," nearest contracts, ",(string count d`changeDates)," change-days):";
-1 "  explained variance per PC (level / slope / curvature): ",-3!d`explainedVar;
{[d;j] -1 "  PC",(string j)," loading: ",-3!(d`loadings)[;j]} [d] each til .cfg.factor`k;
-1 "  (PC1 ~ a flat LEVEL shift, PC2 ~ a SLOPE tilt, PC3 ~ CURVATURE - the textbook term-structure factors.)";
-1 "";

/ (2) the factorRelativeValue template + its own gates.
full:.template.factorRv.run[`commodity`curveHistory!(`CRUDE;clHist)];
stab:full[`validation;`factorStability];
stat:full[`validation;`stationarity];
-1 "factorRelativeValue template gates (template-specific, run BEFORE the universal gates):";
-1 "  factor-stability : cos=",(string stab`cos)," stable=",string stab`stable;
-1 "  residual-stationarity: AR(1)=",(string stat`ar1)," half-life=",(string stat`halfLife)," stationary=",string stat`stationary;
-1 "  spread PnL total : ",string sum (full`pnl)`pnl;
gatesPass:(stab`stable) and stat`stationary;
-1 "";

if[not gatesPass;
    -1 "A template-specific gate FAILED -> the factor-RV thesis is rejected BEFORE the universal";
    -1 "gates run (unstable factors or a non-mean-reverting residual have no tradeable edge)."];

if[gatesPass;
    -1 "Template gates passed -> run the FULL bounded workflow (.workflow.run).";
    runner:{[clHist;from;to] (.template.factorRv.run[`commodity`curveHistory`dateFrom`dateTo!(`CRUDE;clHist;from;to)])`pnl}[clHist];
    .gov.hypoTbl:.gov.__emptyHypotheses[]; .gov.trialTbl:.gov.__emptyTrials[];
    proposal:`hypoId`thesis`edgeSource`instruments`claimedRegimes`capName`runner`axis!(
        `rv_factor_crude;
        "the crude curve reverts toward its k-factor PCA shape (structural)";
        `mispricing;
        `CRUDE;
        `contango`backwardation;
        `factorRelativeValue;
        runner;
        `curveState);
    packet:.workflow.run proposal;
    -1 "  carded=",(string packet`carded),"  trials=",(string .gov.nTrials `rv_factor_crude);
    -1 "  gate verdicts by regime:";
    show select bucket,nObs,netSharpe,dsr,verdict,tradeable from packet`verdicts;
    -1 "  skeptic risk memory: ",$[0<count packet`riskMemory; packet`riskMemory; "(no regime context)"];
    -1 "  ANY tradeable: ",(string packet`anyTradeable),"  decision: ",(string packet`decision)," (human sign-off required: ",(string packet`humanSignOffRequired),")"];

-1 "";
-1 "(R8: a sophisticated shape - factor-structure relative value - added purely as plug-ins (a";
-1 " capability + a template), card-able and judged by the same honest gates. Its main overfitting";
-1 " risk (the choice of k and the lookback) is exactly what the deflation + walk-forward + holdout";
-1 " cascade exists to catch - the gates do their job; the verdict is whatever it honestly is.)";
exit 0;
