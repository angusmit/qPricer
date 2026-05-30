\l core/init.q
/ Cross-check the strategy P&L against the gamma+theta attribution identity. For small
/ spot moves and a fine FDM grid, gammaReconResidual = (totalPnl - txnCostTotal -
/ financingTotal) - (theoreticalGammaPnlTotal + thetaPnlTotal) should be small
/ relative to the option's initial value.
/ Tolerance (documented): residual / optionPrice0 must be within 0.10 with the
/ parameters below. Tighter convergence would need finer FDM grids than are
/ practical in CI.

pathCfg:`spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.10;13;1f%2520f;0f;0f;101);
pathTbl:.strategy.path.fromSynthetic pathCfg;

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `T_RECON;`X;`equityOption;`european;`call;100f;0.5;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;200;200;0f;200f;`linear;1b;1b);

stratCfg:.strategy.defaultConfig `gammaScalp;
stratCfg:@[stratCfg;(`rebalanceMode;`rebalanceInterval;`stepYears;`txnCostRate;`financingRate);:;(`interval;1;1f%2520f;0f;0f)];

bundle:.strategy.runAndSummarize[`gammaScalp;trade;pathTbl;bsModel;fdmCfg;stratCfg];
summary:bundle`summary;
resultTbl:bundle`result;

optionPrice0:resultTbl[0;`optionPrice];
absResidual:abs summary`gammaReconResidual;
relativeResidual:absResidual%optionPrice0;

.testutil.assertTrue[`OK=summary`status;"summary status OK"];
.testutil.assertTrue[not null summary`gammaReconResidual;"residual finite"];
.testutil.assertTrue[relativeResidual<0.10;"residual / optionPrice0 within 0.10 for small-step path"];
.testutil.assertNear[summary`txnCostTotal;0f;1e-12;"zero txnCostRate -> zero txnCostTotal"];
.testutil.assertNear[summary`financingTotal;0f;1e-12;"zero financingRate -> zero financingTotal"];

-1 "PASS test_strategy_gamma_scalp_reconciliation: residual=",string[summary`gammaReconResidual],", optionPrice0=",string[optionPrice0],", relative=",string[relativeResidual];
