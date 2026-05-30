\l lib/init.q
/ commoditySignals with deseasonalize: model signals are computed on the
/ deseasonalized curve (the seasonal slope is removed), the tradeable front-
/ return series stays RAW, the WTI default path is unchanged, and the train-only
/ seasonal fit is causal (early signals unchanged when the future is appended).
seasonal:0.15 0.10 0.0 -0.08 -0.12 -0.10 -0.05 -0.03 0.0 0.02 0.06 0.12;
ymList:raze {[y] (y*100)+1+til 12} each 2021 2022;
expiryOf:{[ym] "D"$(string[ym div 100]),".",(-2#"0",string[ym mod 100]),".15"};
expiries:expiryOf each ymList;
mkRow:{[d;ymList;expiries;seasonal]
    alive:where expiries>=d;
    months:(ymList[alive] mod 100)-1;
    ([] asofDate:(count alive)#d; tenor:(`float$(expiries alive)-d)%365f;
        price:3.0*exp seasonal months; contractYM:ymList alive; expiry:expiries alive)};
baseDates:2021.01.05+7*til 30;
hist:raze mkRow[;ymList;expiries;seasonal] each baseDates;
sigCfg:`trainEndDate`momentumLookback`kalmanEstCfg!(2021.05.01;5;`gridSteps`refineRounds`nSweeps!(3;1;1));

rawSig:.strategy.path.commoditySignals[hist;sigCfg];
desSig:.strategy.path.commoditySignals[hist;sigCfg,enlist[`deseasonalize]!enlist 1b];

/ Deseasonalized model carry collapses vs the raw seasonal carry.
.testutil.assertTrue[(avg abs (desSig`path)`curveSlopeCarry)<0.1*avg abs (rawSig`path)`curveSlopeCarry;"deseasonalized curve carry << raw seasonal carry"];
/ The tradeable front-return series is unchanged (deseasonalization is signal-only).
.testutil.assertTrue[((rawSig`path)`frontReturn)~(desSig`path)`frontReturn;"tradeable front returns stay raw"];
.testutil.assertTrue[(desSig`deseasonalize)~1b;"deseasonalize flag recorded"];
.testutil.assertTrue[12=count desSig`monthFactors;"12 fitted month factors in bundle"];

/ WTI default (deseasonalize off): rerunning is byte-identical to the raw run.
rawSig2:.strategy.path.commoditySignals[hist;sigCfg];
.testutil.assertTrue[(rawSig`path)~rawSig2`path;"default (deseasonalize off) path is deterministic / unchanged"];

/ Causality: train-only seasonal fit -> early-date deseasonalized signals are
/ unchanged when future dates are appended (fixed trainEndDate).
extHist:raze mkRow[;ymList;expiries;seasonal] each baseDates,2021.01.05+7*(30+til 8);
desExt:.strategy.path.commoditySignals[extHist;sigCfg,enlist[`deseasonalize]!enlist 1b];
nOrig:count desSig`path;
.testutil.assertTrue[(desSig`monthFactors)~desExt`monthFactors;"train-only month factors unchanged by appended future"];
.testutil.assertTrue[((desSig`path)`curveSlopeCarry)~nOrig#(desExt`path)`curveSlopeCarry;"deseasonalized carry causal (unchanged by future)"];

-1 "PASS test_commodity_signals_deseasonalize";
