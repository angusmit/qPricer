\l core/init.q
/ ============================================================================
/ carded_gated_research.q - real-data example (NOT a test), DATA-CONDITIONAL.
/ v0.71 connective wiring, made visible: the layers now inform each other (cycle-free).
/   (1) an UNCARDED capability -> .cards.gatedRun REFUSES (document it first; gates NOT run);
/   (2) a CARDED capability -> the gates run via .gov.runFull, and each verdict now carries
/       the regime risk-memory annotation ("resembles <episode>; here's what killed people");
/   (3) .cards.audit[] now dynamically covers the R6 `template` kind too.
/ Skips if no HDB. exit 0;.
/ ============================================================================
hdbPath:.cfg.paths`hdb;
useHdb:0<count @[{[p] key hsym `$p,"/sym"};hdbPath;{[e] ()}];
if[not useHdb;
    -1 "carded_gated_research SKIPPED - no HDB at ",hdbPath," (run scripts/ingest_hdb.q first).";
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
.gov.register `hypoId`thesis`edgeSource`instruments`claimedRegimes`status!(
    `mom_crude_ts;"crude TS momentum earns a positive excess return";`riskPremium;`CRUDE;enlist `backwardation;`research);

-1 "(1) An UNCARDED capability cannot be gated:";
-1 "    gateReady[blackScholesFdm] = ",(string (.cards.gateReady `blackScholesFdm)`ready)," - ",(.cards.gateReady `blackScholesFdm)`reason;
ref:.cards.gatedRun[`mom_crude_ts;`blackScholesFdm;runner;`curveState];
-1 "    .cards.gatedRun verdict = ",(string first ref`verdict)," (gates NOT run, holdout untouched)";
-1 "";

-1 "(2) A CARDED capability (commoditySignalPath) -> gated, with the regime risk-memory annotation:";
-1 "    gateReady[commoditySignalPath] = ",(string (.cards.gateReady `commoditySignalPath)`ready)," - ",(.cards.gateReady `commoditySignalPath)`reason;
verdicts:.cards.gatedRun[`mom_crude_ts;`commoditySignalPath;runner;`curveState];
show select bucket,nObs,netSharpe,verdict,tradeable from verdicts;
-1 "    risk-memory annotation (informational - changed NO verdict):";
{[r] -1 "      ",(string r`bucket),": ",$[0<count r`riskMemory; r`riskMemory; "(no regime context)"]} each verdicts;
-1 "";

-1 "(3) .cards.audit[] now dynamically covers the R6 template kind:";
a:.cards.audit[];
show select kindCovered:count i by kind from a;
-1 "    template kind present in the audit: ",string `template in a`kind;
-1 "";
-1 "(The wiring, cycle-free: cards -> gov (carded gating), gov -> regime (the skeptic), and the";
-1 " audit walks the registry dynamically. You can no longer gate an undocumented capability, the";
-1 " gates surface the relevant historical failure modes, and a new plug-in kind can't hide.)";
exit 0;
