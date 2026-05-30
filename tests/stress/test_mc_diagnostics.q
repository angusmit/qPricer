/ test_mc_diagnostics.q - MC diagnostics stress
\l core/init.q
bsCall:.validation.blackScholesClosedForm[`call;100f;100f;1f;0.05;0f;0.2];

/ European MC convergence
eurPriceFn:{[pathCountVal]
    localCfg:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(pathCountVal;1;42;0b;0b;0.95);
    .montecarlo.priceEuropeanMC[`call;100f;100f;1f;0.05;0f;0.2;localCfg]
 };
eurConv:.convergence.mcConvergenceTable[eurPriceFn;1000 5000 10000;bsCall];
.testutil.assertTrue[3=count eurConv;"3 European convergence rows"];

/ Asian MC convergence
asianTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`averageType`averagingStyle`observationCount!(
    1;`AAPL;`asianOption;`european;`call;100f;1f;1f;`arithmetic;`discrete;50);
asianMkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
asianPriceFn:{[pathCountVal]
    localCfg:enlist[`mcConfig]!enlist `pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(pathCountVal;50;42;0b;0b;0.95);
    .asian.priceAsianOption[asianTrade;asianMkt;localCfg]
 };
asianConv:.convergence.mcConvergenceTable[asianPriceFn;1000 5000 10000;5.8];
.testutil.assertTrue[3=count asianConv;"3 Asian convergence rows"];
.testutil.assertTrue[all (`OK=){x`status} each asianConv;"all Asian rows OK"];

-1 "PASS test_mc_diagnostics: eurRows=",string[count eurConv],", asianRows=",string count asianConv;
