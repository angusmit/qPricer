\l lib/init.q
/ End-to-end walk-forward on a synthetic curve history: produces a per-strategy
/ aggregate distribution + per-split detail + split definitions, with correct
/ shapes. (Synthetic only - no real CSV.)
contracts:202003 202004 202005 202006 202007 202008;
expiries:2020.03.20 2020.04.20 2020.05.20 2020.06.20 2020.07.20 2020.08.20;
mkRows:{[d;contracts;expiries]
    alive:where expiries>=d;
    tens:(`float$(expiries alive)-d)%365f;
    base:55f-2.5f*tens+0.3f*sin 30f*tens;
    ([] asofDate:(count alive)#d; tenor:tens; price:base; contractYM:contracts alive; expiry:expiries alive)};
hist:raze mkRows[;contracts;expiries] each 2020.01.20+til 40;
strategies:`convenienceYieldCarry`timeSeriesMomentum`curveRelativeValue;
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(`T;`WTI;`equityOption;`european;`call;55f;0.25;1f);
sigCfg:`momentumLookback`kalmanEstCfg!(5;`gridSteps`refineRounds`nSweeps!(3;1;1));
splitCfg:`scheme`trainSpan`testSpan`maxSplits!(`expanding;15;7;5);

wf:.strategy.commodityBT.walkForward[strategies;trade;hist;sigCfg;splitCfg];
nSplits:count wf`splits;
.testutil.assertTrue[nSplits>=2;"at least 2 walk-forward splits"];
.testutil.assertTrue[(`splitId`trainStart`trainEnd`testEnd)~cols wf`splits;"split-definition schema"];
.testutil.assertTrue[all (wf`splits)[`trainEnd]<(wf`splits)[`testEnd];"trainEnd precedes testEnd in every split"];
.testutil.assertTrue[(count[strategies]*nSplits)=count wf`detail;"detail has one row per strategy per split"];
agg:wf`aggregate;
.testutil.assertTrue[(`strategyName`nSplits`nTraded`meanSharpe`sharpeStd`splitsPositive`meanAnnReturn`meanMaxDrawdown)~cols agg;"aggregate schema"];
.testutil.assertTrue[(count strategies)=count agg;"one aggregate row per strategy"];
.testutil.assertTrue[all (agg`nSplits)=nSplits;"each strategy aggregated over all splits"];

-1 "PASS test_walk_forward_run: ",(string nSplits)," splits over ",(string count strategies)," strategies";
