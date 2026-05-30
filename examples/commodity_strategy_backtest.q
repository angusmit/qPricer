\l lib/init.q
/ ============================================================================
/ commodity_strategy_backtest.q - real-data example (NOT a test).
/ ----------------------------------------------------------------------------
/ Backtests the v0.53 commodity strategies on the REAL WTI curve history: builds
/ the signal-augmented path (roll-adjusted returns + causal convenience-yield /
/ Kalman-chi / momentum signals, with a fixed train/test split), runs the three
/ headline strategies, and prints the RANKED out-of-sample table + the cross-
/ strategy correlation matrix. Answers: do the calibration/Kalman signals trade?
/ CSVs are user-supplied and gitignored. Non-positive prices are excluded.
/ ============================================================================
crudeDir:"data/barchart/CRUDE";
longTable:.parser.crude.loadAll crudeDir;
allDates:asc distinct longTable`date;
tradeDates:allDates where (allDates>=2020.01.02) and allDates<=2021.12.31;
-1 "Building daily curve history over ",(string count tradeDates)," WTI trading dates ...";
curveHistory:.parser.crude.curveHistory[longTable;tradeDates];

/ Train on 2020 (the regime-rich COVID year), trade out-of-sample on 2021.
trainEnd:2020.12.31;
sigCfg:`rollDaysBeforeExpiry`trainEndDate`momentumLookback`txnCostRate`targetVol`riskFreeRate`kalmanEstCfg!(
    5;trainEnd;20;0.0005;0.15;0.02;`gridSteps`refineRounds`nSweeps!(7;3;3));
-1 "Building signals (Kalman estimated ONCE on train, filtered forward) ...";
sig:.strategy.path.commoditySignals[curveHistory;sigCfg];
p:sig`path;
-1 "Path: ",(string count p)," dates; train ",(string sig`nTrain)," / test ",(string sig`nTest);
-1 "Kalman (train) kappa=",(string (sig`kalmanParams)`kappa)," sigChi=",(string (sig`kalmanParams)`sigChi)," sigXi=",string (sig`kalmanParams)`sigXi;

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `BT;`WTI;`equityOption;`european;`call;60f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:()!();
allStrategies:`convenienceYieldCarry`chiReversion`timeSeriesMomentum`twoTimescale`storageCashCarry`carryMomentumCombo`curveRelativeValue;
suite:.strategy.commodityBT.runSuite[allStrategies;trade;p;bsModel;fdmCfg;()!()];

-1 "";
-1 "Out-of-sample (2021) ranked performance:";
show select strategyName,testSharpe,testAnnualReturn,testAnnualVol,testMaxDrawdown,testHitRate from suite`ranked;
-1 "";
-1 "Cross-strategy out-of-sample P&L correlation:";
-1 "  names : ",-3!suite`correlationNames;
-1 "  matrix: ",-3!suite`correlationMatrix;

/ Calendar-roll harvested yield by regime (a separate, physical curve strategy).
-1 "";
rcc:.strategy.realCurveCalendarRoll[curveHistory;`nearTenor`farTenor`rollTriggerTenor`stepYears!(0.08;0.25;0.16;1f%252f)];
-1 "realCurveCalendarRoll roll yield by regime:";
show rcc`rollByRegime;
-1 "";
-1 "Read: vol-targeted to a common 15% so the ranking compares EDGE not leverage;";
-1 "costs included; signals are causal and the Sharpe is out-of-sample (2021).";
