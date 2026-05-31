\l core/init.q
/ ============================================================================
/ regime_conditional_backtest.q - real-data example (NOT a test), DATA-CONDITIONAL.
/ Research OS R1 (ARCHITECTURE.md 11.1): run an EXISTING commodity strategy through
/ the EXISTING backtest, then break its daily PnL down BY REGIME (curve state) with
/ .regime.breakdown. The strategy's headline is UNCHANGED - the same backtest; only
/ the reporting view is new ("blended 0.8 Sharpe" -> "+X in backwardation / -Y in
/ contango", the tradeable truth). Skips gracefully if the HDB is absent.
/ ============================================================================
hdbPath:.cfg.paths`hdb;
useHdb:0<count @[{[p] key hsym `$p,"/sym"};hdbPath;{[e] ()}];
if[not useHdb;
    -1 "regime_conditional_backtest SKIPPED - no HDB at ",hdbPath," (run scripts/ingest_hdb.q first).";
    exit 0];
.data.hdb.open hdbPath;

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(`RC;`WTI;`equityOption;`european;`call;60f;0.25;1f);
sigCfg:`rollDaysBeforeExpiry`momentumLookback`txnCostRate`targetVol`riskFreeRate`kalmanEstCfg!(
    5;20;0.0005;0.15;0.02;`gridSteps`refineRounds`nSweeps!(7;3;3));
bsModel:.model.createBlackScholesModel[]; fc:()!();

clHist:.data.hdb.curveHistory[`CRUDE;.data.hdb.dates `CRUDE];
sig:.strategy.path.commoditySignals[clHist;sigCfg];
strat:`timeSeriesMomentum;
res:(.strategy.runAndSummarize[strat;trade;sig`path;bsModel;fc;.strategy.defaultConfig strat])`result;

/ the backtest output, unchanged: per-day PnL by date.
ok:select from res where status=`OK;
pnlByDate:select date:stepDate, pnl:stepPnl from ok;
-1 (string strat)," on CRUDE: ",(string count pnlByDate)," days, all-data Sharpe=",
    string (.strategy.commodityBT.__perf[ok`stepPnl;252f;1f])`sharpe;

/ the new view: regime fingerprint per date, then break the SAME PnL down by curve state.
labs:.regime.series[`CRUDE;pnlByDate`date];
br:.regime.breakdown[pnlByDate;labs;`curveState];
-1 "";
-1 "Regime-conditional performance by curve state (the (blended) row == the headline above):";
show select bucket,nDays,sharpe,annualReturn,maxDrawdown,hitRate from br;
-1 "";
-1 "(One blended Sharpe hides whether the edge is structural or one regime - this is the honest view R1 unlocks.)";
exit 0;
