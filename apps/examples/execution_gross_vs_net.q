\l core/init.q
/ ============================================================================
/ execution_gross_vs_net.q - real-data example (NOT a test), DATA-CONDITIONAL.
/ ----------------------------------------------------------------------------
/ The honest deliverable of the v0.60 execution layer: how much of each commodity
/ strategy's GROSS (frictionless) edge survives realistic execution costs. Builds
/ the CRUDE signal path from the HDB once, then runs the strategy suite at a
/ frictionless setting and at several slippage levels, printing GROSS vs NET
/ out-of-sample Sharpe per strategy. Slippage cost = |filled turnover| * refPrice
/ * bps/1e4; because positions are vol-targeted and traded frequently, even a few
/ bps is a material drag - which is exactly the point this layer measures.
/ Skips gracefully if the HDB is absent (build it with scripts/ingest_hdb.q).
/ ============================================================================
hdbPath:.cfg.paths`hdb;
useHdb:0<count @[{[p] key hsym `$p,"/sym"};hdbPath;{[e] ()}];
if[not useHdb;
    -1 "execution_gross_vs_net SKIPPED - no HDB at ",hdbPath," (run scripts/ingest_hdb.q first).";
    exit 0];
.data.hdb.open hdbPath;

strategies:`convenienceYieldCarry`chiReversion`timeSeriesMomentum`twoTimescale`storageCashCarry`carryMomentumCombo`curveRelativeValue;
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(`X;`WTI;`equityOption;`european;`call;60f;0.25;1f);
sigCfg:`rollDaysBeforeExpiry`momentumLookback`txnCostRate`targetVol`riskFreeRate`kalmanEstCfg!(
    5;20;0.0005;0.15;0.02;`gridSteps`refineRounds`nSweeps!(7;3;3));
bsModel:.model.createBlackScholesModel[]; fc:()!();

clHist:.data.hdb.curveHistory[`CRUDE;.data.hdb.dates `CRUDE];
sig:.strategy.path.commoditySignals[clHist;sigCfg];
path:sig`path;
-1 "CRUDE signal path: ",(string count path)," dates ",(string first path`stepDate)," .. ",string last path`stepDate;

/ One strategy's cfg with a slippage-bps execution sub-config; then the cfgByName
/ dict the suite expects (strategy -> cfg). Uses the script globals.
mkExecCfg:{[nm;bps] (.strategy.defaultConfig nm),(enlist `exec)!enlist (enlist `slippageBps)!enlist bps};
mkCfgByName:{[bps] strategies!mkExecCfg[;bps] each strategies};
netSharpeAt:{[bps]
    suite:.strategy.commodityBT.runSuite[strategies;trade;path;bsModel;fc;mkCfgByName bps];
    exec testSharpe by strategyName from suite`ranked};

bpsLevels:0 0.5 1 2 5f;
sharpeByBps:netSharpeAt each bpsLevels;
/ One row per slippage level (bps=0 is GROSS / frictionless), columns = strategies.
report:{[bps;d;strats] (enlist[`slippageBps]!enlist bps),strats!d strats}[;;strategies]'[bpsLevels;sharpeByBps];
-1 "";
-1 "Out-of-sample Sharpe by slippage level (CRUDE, single history); slippageBps=0 is GROSS:";
show report;
-1 "";
gross:first (),(sharpeByBps 0)`timeSeriesMomentum;
net2:first (),(sharpeByBps 3)`timeSeriesMomentum;
-1 "timeSeriesMomentum: gross OOS Sharpe=",(string gross),"  ->  net @2bps=",string net2;
-1 "(The cross-commodity walk-forward verdict of momentum +0.74 is GROSS; this shows the per-bps net decay on CRUDE.)";
exit 0;
