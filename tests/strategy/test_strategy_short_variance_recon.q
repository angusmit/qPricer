\l lib/init.q
/ Reconciliation: for small spot moves and a fine FDM grid, the short-variance
/ gammaReconResidual = (totalPnl - txnCost - financing) - (gammaPnlTotal + thetaPnlTotal)
/ should be small relative to the entry premium. Same tolerance convention as the
/ v0.42 gamma-scalp reconciliation test (documented; valid only for small steps).

pathCfg:`spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.10;13;1f%2520f;0f;0f;101);
pathTbl:.strategy.path.fromSynthetic pathCfg;
pathTbl:update volatility:count[pathTbl]#0.30 from pathTbl;

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `SV_RECON;`X;`equityOption;`european;`call;100f;0.5;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;200;200;0f;200f;`linear;1b;1b);
stratCfg:.strategy.defaultConfig `shortVariance;
stratCfg:@[stratCfg;(`forecastVol;`entryMargin;`rebalanceMode;`rebalanceInterval;`stepYears;`txnCostRate;`financingRate);:;(0.05;0.02;`interval;1;1f%2520f;0f;0f)];

bundle:.strategy.runAndSummarize[`shortVariance;trade;pathTbl;bsModel;fdmCfg;stratCfg];
summary:bundle`summary;
.testutil.assertTrue[summary`gateOpen;"gate open under these params (IV 0.30 > 0.05+0.02)"];
.testutil.assertTrue[`OK=summary`status;"summary status OK"];

absResidual:abs summary`gammaReconResidual;
relResidual:absResidual%summary`premiumCollected;
.testutil.assertTrue[not null summary`gammaReconResidual;"residual finite"];
.testutil.assertTrue[relResidual<0.10;"residual / premiumCollected within 0.10 for small-step path"];
.testutil.assertNear[summary`txnCostTotal;0f;1e-12;"zero txnCostRate -> zero txnCostTotal"];
.testutil.assertNear[summary`financingTotal;0f;1e-12;"zero financingRate -> zero financingTotal"];

.testutil.assertTrue[(summary`theoreticalGammaPnlTotal)<=0f;"short gamma -> negative gamma pnl total"];
.testutil.assertTrue[(summary`thetaPnlTotal)>=0f;"short theta -> non-negative theta pnl total"];

-1 "PASS test_strategy_short_variance_recon: residual=",string[summary`gammaReconResidual],", premium=",string[summary`premiumCollected],", relative=",string[relResidual];
