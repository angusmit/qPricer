\l lib/init.q
/ Margrabe exchange option: payoff max(F1-F2,0) for call. F1=75, F2=70, T=0.25,
/ v1=0.30, v2=0.25, rho=0.8, r=0.05. No strike key (exchange strike is zero).
margrabeParams:`fwd1`fwd2`expiry`vol1`vol2`correlation`riskFreeRate!(75f;70f;0.25;0.30;0.25;0.8;0.05);
mCall:.commodity.spread.margrabePrice[`call;margrabeParams];
mPut:.commodity.spread.margrabePrice[`put;margrabeParams];
.testutil.assertTrue[mCall>0f;"margrabe call positive"];
.testutil.assertTrue[mPut>0f;"margrabe put positive"];

/ Benchmark: Margrabe call MUST equal the Kirk spread call at strike 0 (Kirk
/ reduces exactly to Margrabe when K=0 -> adjStrike=F2, vol2Adj=vol2).
kirkParams0:`fwd1`fwd2`strike`expiry`vol1`vol2`correlation`riskFreeRate!(75f;70f;0f;0.25;0.30;0.25;0.8;0.05);
kirk0:.commodity.spread.kirkPrice[`call;kirkParams0];
.testutil.assertNear[mCall;kirk0;1e-10;"margrabe call == kirk call at strike 0"];

/ Forward parity: call - put = discount * (F1 - F2).
disc:exp neg 0.05*0.25;
.testutil.assertNear[mCall-mPut;disc*75f-70f;1e-10;"margrabe forward parity"];

/ Zero correlation increases spread vol vs high correlation -> higher option value.
lowCorrParams:@[margrabeParams;`correlation;:;0f];
mCallLowCorr:.commodity.spread.margrabePrice[`call;lowCorrParams];
.testutil.assertTrue[mCallLowCorr>mCall;"lower correlation -> higher exchange-option value"];

/ Invalid correlation rejected.
badParams:@[margrabeParams;`correlation;:;1.5];
badRes:@[.commodity.spread.margrabePrice[`call;];badParams;{`ERROR}];
.testutil.assertTrue[badRes~`ERROR;"correlation > 1 rejected"];

-1 "PASS test_spread_option_margrabe: call=",(string mCall),", put=",string mPut;
