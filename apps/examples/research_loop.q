\l core/init.q
/ ============================================================================
/ research_loop.q - real-data example (NOT a test), DATA-CONDITIONAL. Research OS R7.
/ The bounded research loop end-to-end via .workflow.run: a proposal is registered, the
/ capability's carded-ness is checked, the gates run (one-shot holdout) with the skeptic's
/ regime risk-memory annotation and the trial logged, and a HUMAN-ESCALATION PACKET is
/ returned. The HUMAN decides - the loop trades nothing and allocates nothing. The verdict
/ is whatever it honestly is (not tuned). Skips if no HDB. exit 0;.
/ ============================================================================
hdbPath:.cfg.paths`hdb;
useHdb:0<count @[{[p] key hsym `$p,"/sym"};hdbPath;{[e] ()}];
if[not useHdb;
    -1 "research_loop SKIPPED - no HDB at ",hdbPath," (run scripts/ingest_hdb.q first).";
    exit 0];
.data.hdb.open hdbPath;
if[0<count key hsym `$hdbPath,"/regimeEpisodes/"; .regime.library.open hdbPath];

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(`RC;`WTI;`equityOption;`european;`call;60f;0.25;1f);
sigCfg:`rollDaysBeforeExpiry`momentumLookback`txnCostRate`targetVol`riskFreeRate`kalmanEstCfg!(
    5;20;0.0005;0.15;0.02;`gridSteps`refineRounds`nSweeps!(7;3;3));
bsModel:.model.createBlackScholesModel[];
stratCfg:(.strategy.defaultConfig `timeSeriesMomentum),(enlist `exec)!enlist .cfg.gov`exec;
clHist:.data.hdb.curveHistory[`CRUDE;.data.hdb.dates `CRUDE];
ctx:`clHist`sigCfg`bsModel`trade`stratCfg!(clHist;sigCfg;bsModel;trade;stratCfg);
runner:{[ctx;from;to]
    sub:select from ctx[`clHist] where asofDate within (from;to);
    sig:.strategy.path.commoditySignals[sub;ctx`sigCfg];
    res:(.strategy.runAndSummarize[`timeSeriesMomentum;ctx`trade;sig`path;ctx`bsModel;()!();ctx`stratCfg])`result;
    ok:select from res where status=`OK;
    select date:stepDate, pnl:stepPnl from ok}[ctx];

.gov.hypoTbl:.gov.__emptyHypotheses[]; .gov.trialTbl:.gov.__emptyTrials[];

/ the RESEARCHER's proposal (pre-registered economic thesis + claimed edge source).
proposal:`hypoId`thesis`edgeSource`instruments`claimedRegimes`capName`runner`axis!(
    `mom_crude_ts;
    "crude time-series momentum earns a positive excess return (tight-supply trend persistence)";
    `riskPremium;
    `CRUDE;
    enlist `backwardation;
    `commoditySignalPath;
    runner;
    `curveState);

-1 "RESEARCH LOOP for hypothesis ",string proposal`hypoId;
-1 "================================================";
packet:.workflow.run proposal;

-1 "capability   : ",(string packet`capability)," (carded=",(string packet`carded),", gatesRun=",(string packet`gatesRun),")";
-1 "trials logged: ",string .gov.nTrials proposal`hypoId;
-1 "";
if[packet`gatesRun;
    -1 "Gate verdicts by regime (validator + the R3b cascade):";
    show select bucket,nObs,netSharpe,dsr,verdict,tradeable from packet`verdicts;
    -1 "";
    -1 "SKEPTIC risk-memory annotation: ",$[0<count packet`riskMemory; packet`riskMemory; "(no regime context)"]];
-1 "";
-1 "================ HUMAN-ESCALATION PACKET ================";
-1 "anyTradeable (cleared ALL gates in some regime?): ",string packet`anyTradeable;
-1 "decision               : ",string packet`decision;
-1 "humanSignOffRequired   : ",string packet`humanSignOffRequired;
-1 "capital staging        : ",(string packet`stage)," (advances paper->small->scaled ONLY on human sign-off)";
-1 "summary                : ",packet`summary;
-1 "";
-1 "(The loop PROPOSED, gated, challenged, recorded, and PACKAGED - it decided nothing about";
-1 " capital or deployment. The system enforces the discipline: holdout sealed, ledger honest,";
-1 " nothing tradeable without all gates, and the HUMAN holds the go/no-go.)";
exit 0;
