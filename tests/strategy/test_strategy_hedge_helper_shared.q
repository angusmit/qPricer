\l core/init.q
/ Unit-test .strategy.__hedgeInit and .strategy.__hedgeStep on hand-computed inputs.
/ Also re-run a known gamma-scalp scenario to confirm the post-refactor outputs are
/ identical to the v0.42 reference numbers in the PASS line of
/ tests/strategy/test_strategy_gamma_scalp_synthetic.q
/ (steps=11, totalPnl=0.6521844, gammaTot=0.5607279).

initOut:.strategy.__hedgeInit `spot`positionDelta`txnCostRate!(100f;0.5;0.001);
.testutil.assertNear[initOut`hedgePosition;-0.5;1e-12;"hedgePosition = -positionDelta"];
.testutil.assertNear[initOut`hedgeTrade;-0.5;1e-12;"hedgeTrade = hedgePosition at init"];
.testutil.assertNear[initOut`txnCost;(abs -0.5)*100f*0.001;1e-12;"txnCost = |trade|*spot*rate"];
.testutil.assertNear[initOut`cashAdj;(neg -0.5*100f)-initOut`txnCost;1e-12;"cashAdj = -hedgePos*spot - txnCost"];

hedgeStateIn:`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue!(
    49.95; -0.5; 0.5; 1; 100f; 6.5);
stepInputs:`spot`positionValue`positionDelta`stepIndex`stepYears`txnCostRate`financingRate`rebalanceMode`rebalanceInterval`deltaBand!(
    101f; 7.1; 0.55; 1; 1f%252f; 0.001; 0f; `interval; 1; 0.05);
stepOut:.strategy.__hedgeStep[hedgeStateIn;stepInputs];

.testutil.assertNear[stepOut`positionPnl;7.1-6.5;1e-12;"positionPnl = positionValue - prevPositionValue"];
.testutil.assertNear[stepOut`hedgePnl;-0.5*1f;1e-12;"hedgePnl = prevHedgePos*(spot-prevSpot)"];
.testutil.assertNear[stepOut`hedgeTrade;-0.05;1e-12;"hedgeTrade = -0.55 - -0.5 = -0.05; full rebalance to new delta"];
.testutil.assertNear[stepOut`txnCost;(abs stepOut`hedgeTrade)*101f*0.001;1e-12;"txnCost arithmetic"];
.testutil.assertNear[stepOut`stepPnl;(stepOut[`positionPnl]+stepOut`hedgePnl)-stepOut`txnCost;1e-12;"stepPnl with zero financing"];
.testutil.assertNear[stepOut`hedgePosition;-0.55;1e-12;"hedgePosition updated to -positionDelta"];
.testutil.assertNear[stepOut`hedgedDelta;0.55;1e-12;"hedgedDelta updated to positionDelta"];
.testutil.assertTrue[2=stepOut`numRebalances;"numRebalances incremented to 2"];

bandInputs:@[stepInputs;(`rebalanceMode;`deltaBand);:;(`band;0.10)];
bandOut:.strategy.__hedgeStep[hedgeStateIn;bandInputs];
.testutil.assertNear[bandOut`hedgeTrade;0f;1e-12;"band trigger NOT met (|0.55-0.5|=0.05 < 0.10) -> no trade"];
.testutil.assertTrue[1=bandOut`numRebalances;"numRebalances unchanged when no trade"];

pathCfg:`spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.30;11;1f%252f;0.05;0f;42);
pathTbl:.strategy.path.fromSynthetic pathCfg;
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `T_REFAC;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;100;200;0f;300f;`linear;1b;1b);
stratCfg:.strategy.defaultConfig `gammaScalp;
stratCfg:@[stratCfg;(`rebalanceMode;`rebalanceInterval;`stepYears);:;(`interval;1;1f%252f)];
bundle:.strategy.runAndSummarize[`gammaScalp;trade;pathTbl;bsModel;fdmCfg;stratCfg];
summary:bundle`summary;
.testutil.assertNear[summary`totalPnl;0.6521844;1e-6;"gammaScalp totalPnl matches v0.42 reference"];
.testutil.assertNear[summary`theoreticalGammaPnlTotal;0.5607279;1e-6;"gammaScalp gammaTot matches v0.42 reference"];

-1 "PASS test_strategy_hedge_helper_shared: gammaScalpTotalPnl=",string[summary`totalPnl],", initCashAdj=",string[initOut`cashAdj];
