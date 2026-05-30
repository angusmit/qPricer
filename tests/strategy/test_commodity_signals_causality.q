\l lib/init.q
/ Causality / no look-ahead: with a FIXED trainEndDate, every signal at the
/ original dates is UNCHANGED when future dates are appended (the Kalman params
/ come from the fixed train split, the filter is causal, and per-date / trailing
/ signals never see the future). Also checks the train/test split.
contracts:202003 202004 202005 202006 202007;
expiries:2020.03.20 2020.04.20 2020.05.20 2020.06.20 2020.07.20;
mkRows:{[d;contracts;expiries]
    alive:where expiries>=d;
    tens:(`float$(expiries alive)-d)%365f;
    base:55f-2.5f*tens;
    ([] asofDate:(count alive)#d; tenor:tens; price:base; contractYM:contracts alive; expiry:expiries alive)};
baseDates:2020.01.20+til 18;
extDates:2020.01.20+til 26;
baseHist:raze mkRows[;contracts;expiries] each baseDates;
extHist:raze mkRows[;contracts;expiries] each extDates;

sigCfg:`rollDaysBeforeExpiry`trainEndDate`momentumLookback`kalmanEstCfg!(
    5;2020.01.30;3;`gridSteps`refineRounds`nSweeps!(3;1;1));
outBase:.strategy.path.commoditySignals[baseHist;sigCfg];
outExt:.strategy.path.commoditySignals[extHist;sigCfg];
pBase:outBase`path;
nOrig:count pBase;
pExtHead:nOrig#outExt`path;

/ Train/test split is correct against the fixed trainEndDate.
.testutil.assertTrue[(outBase`nTrain)=sum baseDates<=2020.01.30;"train count matches fixed trainEndDate"];
.testutil.assertTrue[(outBase`nTrain)+(outBase`nTest)=count baseDates;"train + test = all dates"];

/ Causality: original-date signals are identical when the future is appended.
.testutil.assertTrue[(pBase`frontReturn)~pExtHead`frontReturn;"frontReturn causal (unchanged by appended future)"];
.testutil.assertTrue[(pBase`momentum)~pExtHead`momentum;"momentum causal"];
.testutil.assertTrue[(pBase`convenienceYield)~pExtHead`convenienceYield;"convenienceYield causal"];
.testutil.assertTrue[(pBase`chi)~pExtHead`chi;"Kalman filtered chi causal (forward filter, fixed-train params)"];
.testutil.assertTrue[(pBase`chiZ)~pExtHead`chiZ;"chiZ causal"];
.testutil.assertTrue[(outBase`volTargetScale)~outExt`volTargetScale;"vol-target scale from fixed train is unchanged"];

-1 "PASS test_commodity_signals_causality: nTrain=",(string outBase`nTrain),", nTest=",string outBase`nTest;
