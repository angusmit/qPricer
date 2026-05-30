\l core/init.q
/ ============================================================================
/ portfolio_allocation.q - real-data example (NOT a test), DATA-CONDITIONAL.
/ Migration step 5 (v0.62): does optimizing the MIX of commodity strategies beat
/ 1/N out-of-sample? Builds each strategy's NET-of-execution daily return series on
/ real WTI (via the HDB), then runs .alloc.compare (causal walk-forward) across the
/ allocation methods. The honest deliverable is the OOS method-comparison table with
/ equalWeight as the baseline. Skips gracefully if the HDB is absent.
/ ============================================================================
hdbPath:.cfg.paths`hdb;
useHdb:0<count @[{[p] key hsym `$p,"/sym"};hdbPath;{[e] ()}];
if[not useHdb;
    -1 "portfolio_allocation SKIPPED - no HDB at ",hdbPath," (run scripts/ingest_hdb.q first).";
    exit 0];
.data.hdb.open hdbPath;

strategies:`convenienceYieldCarry`chiReversion`timeSeriesMomentum`twoTimescale`storageCashCarry`carryMomentumCombo`curveRelativeValue;
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(`PA;`WTI;`equityOption;`european;`call;60f;0.25;1f);
sigCfg:`rollDaysBeforeExpiry`momentumLookback`txnCostRate`targetVol`riskFreeRate`kalmanEstCfg!(
    5;20;0.0005;0.15;0.02;`gridSteps`refineRounds`nSweeps!(7;3;3));
bsModel:.model.createBlackScholesModel[]; fc:()!();

clHist:.data.hdb.curveHistory[`CRUDE;.data.hdb.dates `CRUDE];
sig:.strategy.path.commoditySignals[clHist;sigCfg];
path:sig`path;

/ each strategy's full NET-of-execution daily return series (aligned, same path).
netReturn:{[trade;path;m;fc;nm]
    res:(.strategy.runAndSummarize[nm;trade;path;m;fc;.strategy.defaultConfig nm])`result;
    exec stepPnl from res where status=`OK}[trade;path;bsModel;fc];
panel:netReturn each strategies;
-1 "Strategy return panel: ",(string count strategies)," strategies x ",(string count first panel)," dates";

splitCfg:`scheme`trainSpan`testSpan`maxSplits!(`rolling;252;63;10);
methods:`equalWeight`inverseVol`minVariance`riskParity`maxSharpe`meanVariance;
cmp:.alloc.compare[panel;methods;.alloc.defaultConfig[];splitCfg];
-1 "";
-1 "OOS portfolio comparison (",(string first cmp`nSplits)," walk-forward splits; equalWeight is the baseline):";
show select method,oosSharpe,oosAnnReturn,oosMaxDrawdown,avgTurnover from cmp;
-1 "";
ewSharpe:first exec oosSharpe from cmp where method=`equalWeight;
rpSharpe:first exec oosSharpe from cmp where method=`riskParity;
-1 "Does optimization beat 1/N OOS?  equalWeight Sharpe=",(string ewSharpe),"  riskParity=",(string rpSharpe);
-1 "(riskParity/inverseVol use only the covariance; maxSharpe/meanVariance use return forecasts and tend to overfit OOS - the honest finding.)";
exit 0;
