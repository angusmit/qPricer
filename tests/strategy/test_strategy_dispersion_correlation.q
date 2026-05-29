\l lib/init.q
/ Self-correlation identity check: with all weights equal and a known correlation
/ matrix, the implied avg correlation backed out from the indexVol formula should
/ recover the input average pairwise correlation when indexVol is computed from
/ the analytic formula sigma_idx^2 = w'Sigma w.
namesL:`A`B;
weightsL:0.5 0.5;
volsL:0.20 0.30;
inputAvgCorr:0.4;
corrMatL:(1 0.4f;0.4 1f);
/ Construct indexVol = sqrt(w' Sigma w)
/ Sigma_ij = rho_ij * sig_i * sig_j
sig00:volsL[0]*volsL[0];
sig11:volsL[1]*volsL[1];
sig01:volsL[0]*volsL[1]*0.4;
indexVarAnalytic:((weightsL[0]*weightsL[0])*sig00)+((weightsL[1]*weightsL[1])*sig11)+2f*weightsL[0]*weightsL[1]*sig01;
indexVolAnalytic:sqrt indexVarAnalytic;
pathCfg:`names`weights`spot0`vols`correlationMatrix`steps`stepYears`riskFreeRate`dividendYields`seed!(
    namesL;weightsL;100 100f;volsL;corrMatL;6;1f%252f;0.02;0 0f;42);
bundle:.strategy.path.fromCorrelated pathCfg;
pathTbl:bundle`indexPath;
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `DS_C;`IDX;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;80;150;0f;300f;`linear;1b;1b);
stratCfg:.strategy.defaultConfig `dispersion;
stratCfg:@[stratCfg;(`pathBundle;`indexVol;`stepYears);:;(bundle;indexVolAnalytic;1f%252f)];
runBundle:.strategy.runAndSummarize[`dispersion;trade;pathTbl;bsModel;fdmCfg;stratCfg];
summary:runBundle`summary;
recoveredCorr:summary`impliedCorrAtEntry;
.testutil.assertTrue[1e-10>abs recoveredCorr-inputAvgCorr;"implied corr backed out matches input avg corr"];
-1 "PASS test_strategy_dispersion_correlation: impliedRecovered=",string[recoveredCorr],", inputAvgCorr=",string[inputAvgCorr];
