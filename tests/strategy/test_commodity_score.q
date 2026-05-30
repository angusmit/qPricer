\l lib/init.q
/ Scoring & ranking: the per-strategy time-series metrics match an INDEPENDENT
/ recomputation, the ranking is a stable descending sort by testSharpe, and the
/ correlation matrix is symmetric with a unit diagonal.
n:12;
dates:2020.01.01+til n;
cy:0.10 0.08 0.06 0.04 0.02 -0.02 -0.04 -0.06 -0.08 -0.06 -0.04 0.02;
chiZ:-2.0 -1.5 -0.1 0.2 1.6 2.0 0.5 -1.2 -0.3 0.4 1.1 -1.5;
ret:0.01 -0.005 0.02 0.01 -0.01 0.015 -0.02 0.005 0.008 -0.012 0.004 0.01;
isT:n#0b; isT[til 6]:1b;
p:flip `stepIndex`stepDate`spot`volatility`riskFreeRate`dividendYield`marketPrice`status`frontReturn`frontPrice`volTargetScale`convenienceYield`curveSlopeCarry`chiZ`isTrain!(
    til n; dates; 60f+sums ret; n#0Nf; n#0.02; n#0f; 60f+sums ret; n#`OK; ret; 60f+sums ret; n#1f; cy; cy; chiZ; isT);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(`T;`WTI;`equityOption;`european;`call;60f;0.25;1f);
m:.model.createBlackScholesModel[]; fc:()!();

suite:.strategy.commodityBT.runSuite[`convenienceYieldCarry`chiReversion;trade;p;m;fc;()!()];
ranked:suite`ranked;
.testutil.assertTrue[2=count ranked;"both strategies scored"];

/ Ranking is a stable descending sort by testSharpe.
.testutil.assertTrue[(ranked`testSharpe)~desc ranked`testSharpe;"ranked descending by testSharpe"];

/ Independent recomputation of carry's out-of-sample Sharpe from its stepPnl.
carryRes:(first suite[`results] where (suite[`results][;`name])=`convenienceYieldCarry)`summary;
carryRun:.strategy.run[`convenienceYieldCarry;trade;p;m;fc;.strategy.defaultConfig `convenienceYieldCarry];
testStepPnl:exec stepPnl from carryRun where status=`OK, not isTrain;
dailyRet:testStepPnl%1f;
indepSharpe:$[0f<dev dailyRet;(252f*avg dailyRet)%(dev dailyRet)*sqrt 252f;0Nf];
.testutil.assertNear[carryRes`testSharpe;indepSharpe;1e-9;"reported testSharpe == independent recomputation"];

/ Correlation matrix: symmetric, unit diagonal.
cm:suite`correlationMatrix;
.testutil.assertNear[cm[0;0];1f;1e-12;"unit diagonal [0;0]"];
.testutil.assertNear[cm[1;1];1f;1e-12;"unit diagonal [1;1]"];
.testutil.assertNear[cm[0;1];cm[1;0];1e-12;"correlation matrix symmetric"];

-1 "PASS test_commodity_score";
