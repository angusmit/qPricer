\l lib/init.q
/ ============================================================================
/ cross_commodity_backtest.q - real-data example (NOT a test), DATA-CONDITIONAL.
/ ----------------------------------------------------------------------------
/ Validates the 8-strategy suite across TWO commodities: an extended WTI walk-
/ forward (all liquid dates) and, IF Henry Hub gas data is present, NG with
/ DESEASONALIZED signals. Aggregates into a cross-commodity robustness verdict
/ (mean OOS Sharpe + dispersion + fraction of (commodity x split) cells positive)
/ and contrasts raw vs deseasonalized gas carry. CSVs are user-supplied/gitignored.
/ ============================================================================
strategies:`convenienceYieldCarry`chiReversion`timeSeriesMomentum`twoTimescale`storageCashCarry`carryMomentumCombo`curveRelativeValue;
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(`XC;`COMMODITY;`equityOption;`european;`call;60f;0.25;1f);
sigCfgBase:`rollDaysBeforeExpiry`momentumLookback`txnCostRate`targetVol`riskFreeRate`kalmanEstCfg!(
    5;20;0.0005;0.15;0.02;`gridSteps`refineRounds`nSweeps!(7;3;3));
splitCfg:`scheme`trainSpan`testSpan`maxSplits!(`rolling;252;63;10);
bsModel:.model.createBlackScholesModel[]; fc:()!();

/ --- extended WTI walk-forward over all liquid dates ---
clLong:.parser.futures.loadAll["data/barchart/CRUDE";`CRUDE];
clDates:asc distinct clLong`date;
clHist:.parser.crude.curveHistory[clLong;clDates];
-1 "WTI: ",(string count clDates)," liquid dates ",(string first clDates)," .. ",string last clDates;
wfCL:.strategy.commodityBT.walkForward[strategies;trade;clHist;sigCfgBase;splitCfg];
-1 "WTI extended walk-forward (",(string count wfCL`splits)," rolling 252/63 splits):";
show `meanSharpe xdesc wfCL`aggregate;

/ --- NG (deseasonalized), data-conditional ---
gasFiles:@[{key hsym `$x};"data/barchart/GAS";{[e] ()}];
if[0=count gasFiles;
    -1 "";-1 "NG cross-commodity leg SKIPPED - no data/barchart/GAS (data-conditional).";
    exit 0];
gasLong:.parser.futures.loadAll["data/barchart/GAS";`GAS];
gasDates:asc distinct gasLong`date;
gasHist:.parser.futures.curveHistory[gasLong;gasDates];
-1 "";-1 "NG: ",(string count gasDates)," liquid dates ",(string first gasDates)," .. ",string last gasDates;
wfNG:.strategy.commodityBT.walkForward[strategies;trade;gasHist;sigCfgBase,enlist[`deseasonalize]!enlist 1b;splitCfg];
-1 "NG (deseasonalized) walk-forward (",(string count wfNG`splits)," splits):";
show `meanSharpe xdesc wfNG`aggregate;

/ --- cross-commodity robustness verdict ---
verdict:.strategy.commodityBT.crossCommodity[`CL`NG!(wfCL`detail;wfNG`detail)];
-1 "";-1 "CROSS-COMMODITY robustness verdict (CL + NG, all splits):";
show verdict;

/ --- raw vs deseasonalized gas CARRY contrast (one mid-history split) ---
midTrain:gasDates (count gasDates) div 2;
ngCfgRaw:sigCfgBase,`trainEndDate`deseasonalize!(midTrain;0b);
ngCfgDes:sigCfgBase,`trainEndDate`deseasonalize!(midTrain;1b);
carrySharpe:{[hist;cfg;trade;m;fc]
    sig:.strategy.path.commoditySignals[hist;cfg];
    s:.strategy.commodityBT.runSuite[enlist `convenienceYieldCarry;trade;sig`path;m;fc;()!()];
    first exec testSharpe from s`ranked}[gasHist;;trade;bsModel;fc];
-1 "";-1 "Gas carry OOS Sharpe - raw=",(string carrySharpe ngCfgRaw)," vs deseasonalized=",string carrySharpe ngCfgDes;
exit 0;
