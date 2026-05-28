\l lib/init.q
/ Kirk spread call: F1=75, F2=70, K=3, T=0.25, v1=0.30, v2=0.25, rho=0.8, r=0.05
spreadParams:`fwd1`fwd2`strike`expiry`vol1`vol2`correlation`riskFreeRate!(75f;70f;3f;0.25;0.30;0.25;0.8;0.05);
kirkCall:.commodity.spread.kirkPrice[`call;spreadParams];
kirkPut:.commodity.spread.kirkPrice[`put;spreadParams];
.testutil.assertTrue[kirkCall>0f;"spread call positive"];
.testutil.assertTrue[kirkPut>0f;"spread put positive"];

/ Spread payoff
payoff:.commodity.spread.spreadPayoff[75f;70f;3f;`call];
.testutil.assertNear[payoff;2f;0.01;"payoff = max(75-70-3,0) = 2"];

/ Invalid correlation
badParams:`fwd1`fwd2`strike`expiry`vol1`vol2`correlation`riskFreeRate!(75f;70f;3f;0.25;0.30;0.25;1.5;0.05);
badCorr:@[.commodity.spread.kirkPrice[`call;];badParams;{`ERROR}];
.testutil.assertTrue[badCorr~`ERROR;"correlation > 1 fails"];

-1 "PASS test_spread_option: call=",string[kirkCall],", put=",string kirkPut;
