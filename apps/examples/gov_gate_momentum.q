\l core/init.q
/ ============================================================================
/ gov_gate_momentum.q - real-data example (NOT a test), DATA-CONDITIONAL.
/ Research OS R3 + R3b (ARCHITECTURE.md 11.4): run R1's regime-conditional momentum
/ finding through the FULL 5-gate cascade (.gov.runFull) and report the HONEST,
/ FAIL-SAFE verdict.
/ -
/ R1 surfaced: crude TS momentum's blended Sharpe ~0.31 decomposes into backwardation,
/ contango, flat. R3 added thesis -> cost -> deflated Sharpe -> walk-forward; R3b adds
/ the SEALED HOLDOUT (Gate 4, one-shot) and a FAIL-SAFE verdict (a gate FAILURE is never
/ tradeable). The cascade:
/   (1) register the hypothesis - claimedRegimes = `backwardation only (the original
/       tight-supply trend story), so contango/flat fire the POST-HOC flag;
/   (2) runner[from;to] wraps the EXISTING momentum backtest restricted to a date range
/       (gov controls which dates the strategy EVER sees) and returns NET daily PnL under
/       a realistic .exec cost;
/   (3) .gov.runFull restricts to TRAIN+VALIDATE, runs gates 0-3 per regime bucket, and
/       ONLY if a bucket clears 0-3 does the hypothesis EARN the one-shot holdout look;
/   (4) print, per bucket, the verdict + failed gate + the `tradeable` flag.
/ EXPECTED (honest, not tuned): the slices that fail cost/deflation NEVER reach the
/ holdout (the seal protects it from findings that have not earned a look) - verdict
/ `reject`/`research`, tradeable=0b. The deflation-failed contango/flat slices are
/ `research` (a real regime question, failed evidence), NOT a tradeable regimeConditional.
/ Skips gracefully if the HDB is absent. (.gov.open/.gov.flush persist the ledger to the
/ HDB; this example keeps an in-memory ledger so N reflects just this run.)
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
bsModel:.model.createBlackScholesModel[];
strat:`timeSeriesMomentum;
stratCfg:(.strategy.defaultConfig strat),(enlist `exec)!enlist .cfg.gov`exec;
clHist:.data.hdb.curveHistory[`CRUDE;.data.hdb.dates `CRUDE];

/ runner[from;to]: the EXISTING momentum backtest restricted to [from;to] -> NET per-day
/ PnL. Wraps the unchanged engine; gov decides which dates it sees by slicing clHist.
ctx:`clHist`sigCfg`bsModel`trade`strat`stratCfg!(clHist;sigCfg;bsModel;trade;strat;stratCfg);
runnerBase:{[ctx;from;to]
    sub:select from ctx[`clHist] where asofDate within (from;to);
    sig:.strategy.path.commoditySignals[sub;ctx`sigCfg];
    res:(.strategy.runAndSummarize[ctx`strat;ctx`trade;sig`path;ctx`bsModel;()!();ctx`stratCfg])`result;
    ok:select from res where status=`OK;
    select date:stepDate, pnl:stepPnl from ok };
runner:runnerBase[ctx];

/ (1) pre-register; claimedRegimes deliberately EXCLUDES contango/flat.
.gov.hypoTbl:.gov.__emptyHypotheses[]; .gov.trialTbl:.gov.__emptyTrials[];
.gov.register `hypoId`thesis`edgeSource`instruments`claimedRegimes`status!(
    `mom_crude_ts;
    "crude time-series momentum earns a positive excess return (tight-supply trend persistence)";
    `riskPremium;
    `CRUDE;
    enlist `backwardation;
    `research);

tvR:.gov.zone.range[`CRUDE;`trainValidate];
hoR:.gov.zone.range[`CRUDE;`holdout];
-1 "Zones (CRUDE): trainValidate ",(string tvR 0)," .. ",(string tvR 1),"   |   SEALED holdout ",(string hoR 0)," .. ",string hoR 1;

/ (3) the FULL fail-safe cascade (gates 0-3 on trainValidate; holdout only if earned).
verdicts:.gov.runFull[`mom_crude_ts;runner;`curveState];

-1 "";
-1 "Full 5-gate cascade by curve state (claimedRegimes = `backwardation):";
show select bucket,nObs,netSharpe,dsr,failedGate,verdict,tradeable,postHoc from verdicts;
-1 "";
-1 "Holdout earned? ",$[any verdicts`tradeable; "yes - a bucket cleared gates 0-3"; "NO - no bucket cleared gates 0-3, so the sealed holdout was never read (the seal held)."];
-1 "holdoutUsedAt on the hypothesis: ",string (.gov.hypo `mom_crude_ts)`holdoutUsedAt;
-1 "";
-1 "Reasons:";
{[v] -1 "  ",(string v`bucket),": ",v`reason} each verdicts;
-1 "";
-1 "(R1 surfaced the structure; R3 gates it on cost + deflation + walk-forward; R3b seals the";
-1 " holdout and makes verdicts fail-safe - a deflation-FAILED slice is `research`, never a";
-1 " tradeable `regimeConditional, and tradeable=0b across the board here.)";
exit 0;
