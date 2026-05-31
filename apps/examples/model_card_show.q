\l core/init.q
/ ============================================================================
/ model_card_show.q - real-data example (NOT a test), DATA-CONDITIONAL. Research OS R5.
/ The synthesis payoff: ONE command shows everything known about a capability - identity
/ + the R2 contract it satisfies + assumptions + named edge source + LIVE validation
/ status DERIVED from the gov ledger (honest, not asserted) + linked regime risk memory
/ + regime applicability - then runs .cards.audit[]. exit 0;.
/ -
/ To make the validation line illustrative rather than `ungated`, it first runs the
/ momentum backtest through .gov.run (gates 0-3) so the ledger holds a REAL verdict for
/ the signal's hypothesis, then shows the card deriving that status. Skips if no HDB.
/ ============================================================================
hdbPath:.cfg.paths`hdb;
useHdb:0<count @[{[p] key hsym `$p,"/sym"};hdbPath;{[e] ()}];
if[not useHdb;
    -1 "model_card_show SKIPPED - no HDB at ",hdbPath," (run scripts/ingest_hdb.q first).";
    exit 0];
.data.hdb.open hdbPath;
if[0<count key hsym `$hdbPath,"/regimeEpisodes/"; .regime.library.open hdbPath];

cap:`commoditySignalPath;

/ --- illustrative: log a REAL governance verdict for this signal's hypothesis ---
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(`RC;`WTI;`equityOption;`european;`call;60f;0.25;1f);
sigCfg:`rollDaysBeforeExpiry`momentumLookback`txnCostRate`targetVol`riskFreeRate`kalmanEstCfg!(
    5;20;0.0005;0.15;0.02;`gridSteps`refineRounds`nSweeps!(7;3;3));
bsModel:.model.createBlackScholesModel[];
stratCfg:(.strategy.defaultConfig `timeSeriesMomentum),(enlist `exec)!enlist .cfg.gov`exec;
clHist:.data.hdb.curveHistory[`CRUDE;.data.hdb.dates `CRUDE];
sig:.strategy.path.commoditySignals[clHist;sigCfg];
res:(.strategy.runAndSummarize[`timeSeriesMomentum;trade;sig`path;bsModel;()!();stratCfg])`result;
ok:select from res where status=`OK;
pnlByDate:select date:stepDate, pnl:stepPnl from ok;
labs:.regime.series[`CRUDE;pnlByDate`date];
.gov.hypoTbl:.gov.__emptyHypotheses[]; .gov.trialTbl:.gov.__emptyTrials[];
.gov.register `hypoId`thesis`edgeSource`instruments`claimedRegimes`status!(
    `mom_crude_ts;"crude TS momentum earns a positive excess return";`riskPremium;`CRUDE;enlist `backwardation;`research);
.gov.run[`mom_crude_ts;pnlByDate;labs;`curveState];

/ --- the card synthesis ---
card:.cards.get cap;
kind:card`capabilityKind;
con:.registry.contract kind;
-1 "================ MODEL CARD: ",(string cap)," ================";
-1 "kind/version : ",(string kind)," / ",(string card`version),"   owner: ",(card`owner),"   asOf: ",string card`asOf;
-1 "intended use : ",card`intendedUse;
-1 "assumptions  : ",card`assumptions;
-1 "edge source  : ",string card`edgeSource;
-1 "regimes      : ",card`regimeApplicability;
-1 "";
-1 "R2 CONTRACT it satisfies (",(string kind)," v",(string con`version),"):";
-1 "  required in : ",(", " sv string key con`requiredIn);
-1 "  required out: ",(", " sv string key con`requiredOut);
-1 "  conforms    : ",string .registry.conforms[kind;cap];
-1 "";
-1 "VALIDATION STATUS (DERIVED from the gov ledger - not asserted):";
vstat:.cards.validationStatus cap;
-1 "  hypoId=",(string vstat`hypoId),"  status=",(string vstat`status),"  N=",(string vstat`nTrials);
-1 "  headlineNetSharpe=",(string vstat`headlineNetSharpe),"  headlineDSR=",(string vstat`headlineDsr),"  deflationPass=",string vstat`deflationPass;
-1 "  holdoutUsed=",(string vstat`holdoutUsed),"  TRADEABLE=",string vstat`tradeable;
-1 "";
-1 "LINKED RISK MEMORY (from the regime library, key=",(string card`riskMemoryKey),"):";
rm:.cards.riskMemory cap;
$[0=count rm; -1 "  (none / library not built)";
    {[r] -1 "  ",(string r`label),": ",r`riskMemory} each rm];
-1 "";
-1 "CARD AUDIT (.cards.audit[]) - coverage across all registered capabilities:";
a:.cards.audit[];
-1 "  ",(string sum a`pass)," documented & conforming / ",(string count a)," registered; ",(string sum not a`pass)," awaiting cards:";
-1 "  undocumented: ",", " sv string exec name from a where not pass;
-1 "";
-1 "(One command: what it IS (R2 contract) + assumes + claims (edge) + where it's valid (R1/R4 risk";
-1 " memory) + what the gates ACTUALLY granted (R3/R3b ledger). The validation line is DERIVED - here";
-1 " momentum is logged but not promoted, so it reads inResearch / not tradeable, never 'validated'.)";
exit 0;
