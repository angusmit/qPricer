\l lib/init.q
/ Generate price with known vol, then recover it
trueVol:0.30;
modelPx:.commodity.black76.price[`call;75f;75f;0.25;trueVol;0.05];
solverCfg:`loVol`hiVol`tol!(0.01;1.0;1e-8);
ivResult:.commodity.black76.impliedVol[`call;modelPx;75f;75f;0.25;0.05;solverCfg];
.testutil.assertTrue[ivResult[`status]=`OK;"IV converged"];
.testutil.assertNear[ivResult`impliedVol;trueVol;0.001;"recovered vol"];

/ Invalid market price
badIv:@[{.commodity.black76.impliedVol[`call;x;75f;75f;0.25;0.05;`loVol`hiVol`tol!(0.01;1.0;1e-8)]};-1f;{`ERROR}];
.testutil.assertTrue[badIv~`ERROR;"negative price fails"];

-1 "PASS test_black76_implied_vol: recoveredVol=",string ivResult`impliedVol;
