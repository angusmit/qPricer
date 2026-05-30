\l core/init.q
/ Known-answer + constraint tests for the .alloc allocator (synthetic, deterministic).
/ Each method matches its closed form / definition; constraints bind; weights are valid.
cfgU:.alloc.defaultConfig[],`longOnly`fullyInvested!(0b;1b);   / unconstrained, fully-invested
r:(1.0 -0.5 0.6 0.2 -0.3 0.4 0.1 -0.2;
   0.2 0.1 -0.4 0.3 0.1 -0.2 0.05 0.15;
   -0.3 0.4 0.2 -0.1 0.25 -0.15 0.1 0.05);
n:count r;
cmat:.alloc.covariance[r;.alloc.defaultConfig[]];
sig:sqrt .alloc.__diag cmat;
mu:avg each r;

.testutil.assertTrue[(.alloc.weights[r;`equalWeight;.alloc.defaultConfig[]])~n#1f%n;"equalWeight == 1/N exactly"];
.testutil.assertTrue[1e-12>max abs (.alloc.weights[r;`inverseVol;cfgU])-((reciprocal sig)%sum reciprocal sig);"inverseVol == 1/sigma normalized"];
mv:(inv cmat) mmu n#1f;
.testutil.assertTrue[1e-10>max abs (.alloc.weights[r;`minVariance;cfgU])-(mv%sum mv);"minVariance == Sigma^-1 1 / (1' Sigma^-1 1)"];
ms:(inv cmat) mmu mu;
.testutil.assertTrue[1e-10>max abs (.alloc.weights[r;`maxSharpe;cfgU])-(ms%sum ms);"maxSharpe == Sigma^-1 mu normalized"];
cfgMV:.alloc.defaultConfig[],`longOnly`fullyInvested!(0b;0b); lam:cfgMV`riskAversion;
.testutil.assertTrue[1e-10>max abs (.alloc.weights[r;`meanVariance;cfgMV])-((1f%lam)*(inv cmat) mmu mu);"meanVariance == (1/lambda) Sigma^-1 mu unconstrained"];

/ riskParity: equal risk contributions, long-only, sum 1.
wrp:.alloc.weights[r;`riskParity;.alloc.defaultConfig[]]; rc:wrp*cmat mmu wrp;
.testutil.assertTrue[1e-6>max[rc]-min rc;"riskParity: equal risk contributions w_i*(Sigma w)_i"];
.testutil.assertTrue[(1e-12>abs 1f-sum wrp) and all wrp>=0f;"riskParity: long-only, sums to 1"];

/ every method under the DEFAULT (long-only, fully-invested) gives valid weights.
methAll:`equalWeight`inverseVol`minVariance`maxSharpe`riskParity`meanVariance;
dfw:.alloc.weights[r;;.alloc.defaultConfig[]] each methAll;
.testutil.assertTrue[all {(1e-9>abs 1f-sum x) and all x>=neg 1e-12} each dfw;"all methods: long-only + sum-1 under default config"];

/ weightCap binds; turnover penalty reduces churn.
wc:.alloc.weights[r;`minVariance;.alloc.defaultConfig[],(enlist `weightCap)!enlist 0.5];
.testutil.assertTrue[all wc<=0.5+1e-9;"weightCap caps individual weights"];
.testutil.assertTrue[1e-9>abs 1f-sum wc;"capped weights still sum to 1"];
prevW:n#1f%n; w0:.alloc.weights[r;`maxSharpe;.alloc.defaultConfig[]];
wk:.alloc.weights[r;`maxSharpe;.alloc.defaultConfig[],`turnoverPenalty`prevWeights!(3f;prevW)];
.testutil.assertTrue[(sum abs wk-prevW)<sum abs w0-prevW;"turnoverPenalty reduces weight churn"];

-1 "PASS test_alloc_weights: methods match closed forms; constraints bind";
