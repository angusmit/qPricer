\l core/init.q
/ ============================================================================
/ gov_gate_momentum.q - real-data example (NOT a test), DATA-CONDITIONAL.
/ Research OS R3 (ARCHITECTURE.md 11.4): run R1's regime-conditional momentum
/ finding through the R3 GATE CASCADE and report the HONEST verdict.
/ -
/ R1 surfaced: crude TS momentum's blended Sharpe ~0.31 decomposes into
/ backwardation -0.08 (~1460 days), contango +0.89 (249), flat +2.56 (157).
/ The +2.56 is TEMPTING and statistically FRAGILE - tiny samples, and slicing into
/ three regime buckets is effectively three trials. R3 refuses to be fooled:
/   (1) register the hypothesis - claimedRegimes = `backwardation only (the original
/       tight-supply trend story), so contango/flat fire the POST-HOC flag;
/   (2) run the backtest ONCE with a REALISTIC .exec cost config (.cfg.gov.exec),
/       take the NET daily PnL, regime-break it (.regime.breakdown's grouping);
/   (3) .gov.run each bucket -> logs a trial per slice (N grows to 3, the honest
/       multiple-testing denominator) THEN evaluates the gates;
/   (4) print the per-bucket verdict + DSR.
/ EXPECTED: the contango / flat buckets FAIL Gate 2 (deflated Sharpe) - small n plus
/ the trial-count penalty drive DSR below 0.95 - and are flagged post-hoc. Thresholds
/ are NOT tuned to force any outcome; the gates report exactly what they compute.
/ Skips gracefully if the HDB is absent. (.gov.open / .gov.flush persist the ledger
/ to the HDB; this example keeps an in-memory ledger so N reflects just this run.)
/ ============================================================================
hdbPath:.cfg.paths`hdb;
useHdb:0<count @[{[p] key hsym `$p,"/sym"};hdbPath;{[e] ()}];
if[not useHdb;
    -1 "gov_gate_momentum SKIPPED - no HDB at ",hdbPath," (run scripts/ingest_hdb.q first).";
    exit 0];
.data.hdb.open hdbPath;

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(`RC;`WTI;`equityOption;`european;`call;60f;0.25;1f);
sigCfg:`rollDaysBeforeExpiry`momentumLookback`txnCostRate`targetVol`riskFreeRate`kalmanEstCfg!(
    5;20;0.0005;0.15;0.02;`gridSteps`refineRounds`nSweeps!(7;3;3));
bsModel:.model.createBlackScholesModel[]; fc:()!();
strat:`timeSeriesMomentum;

clHist:.data.hdb.curveHistory[`CRUDE;.data.hdb.dates `CRUDE];
sig:.strategy.path.commoditySignals[clHist;sigCfg];

/ (2) run ONCE with a REALISTIC execution cost config -> NET daily PnL.
stratCfg:(.strategy.defaultConfig strat),(enlist `exec)!enlist .cfg.gov`exec;
res:(.strategy.runAndSummarize[strat;trade;sig`path;bsModel;fc;stratCfg])`result;
ok:select from res where status=`OK;
pnlByDate:select date:stepDate, pnl:stepPnl from ok;
blended:(.strategy.commodityBT.__perf[ok`stepPnl;252f;1f])`sharpe;
-1 (string strat)," on CRUDE (net of realistic execution cost): ",(string count pnlByDate),
    " days, blended Sharpe=",.gov.__fmt blended;

/ regime labels per date.
labs:.regime.series[`CRUDE;pnlByDate`date];

/ (1) pre-register the hypothesis. claimedRegimes deliberately EXCLUDES contango/flat.
.gov.hypoTbl:.gov.__emptyHypotheses[]; .gov.trialTbl:.gov.__emptyTrials[];
.gov.register `hypoId`thesis`edgeSource`instruments`claimedRegimes`status!(
    `mom_crude_ts;
    "crude time-series momentum earns a positive excess return (tight-supply trend persistence)";
    `riskPremium;
    `CRUDE;
    enlist `backwardation;
    `research);

/ (3) gate every regime bucket (logs a trial per slice -> N grows, THEN evaluates).
verdicts:.gov.run[`mom_crude_ts;pnlByDate;labs;`curveState];

-1 "";
-1 "Gate cascade by curve state (claimedRegimes = `backwardation; N = ",(string .gov.nTrials `mom_crude_ts)," trials in the ledger):";
show select bucket,nObs,netSharpe,dsr,postHoc,verdict,failedGate from verdicts;
-1 "";
-1 "Reasons:";
{[v] -1 "  ",(string v`bucket),": ",v`reason} each verdicts;
-1 "";
-1 "(R1 surfaced the +2.56 flat / +0.89 contango structure; R3 refuses to be fooled by it -";
-1 " small n + the 3-trial deflation penalty + the post-hoc flag sink the tempting slices.)";
exit 0;
