\l core/init.q
/ Black-76: F=75, K=75, T=0.25, vol=0.30, r=0.05
callPx:.commodity.black76.price[`call;75f;75f;0.25;0.30;0.05];
putPx:.commodity.black76.price[`put;75f;75f;0.25;0.30;0.05];
.testutil.assertTrue[callPx>0f;"call positive"];
.testutil.assertTrue[putPx>0f;"put positive"];

/ Put-call parity: call - put = exp(-rT) * (F - K) = 0 (ATM)
discFactor:exp neg 0.05*0.25;
parity:(callPx-putPx)-discFactor*75f-75f;
.testutil.assertNear[parity;0f;0.001;"put-call parity ATM"];

/ OTM call
otmCall:.commodity.black76.price[`call;75f;85f;0.25;0.30;0.05];
.testutil.assertTrue[otmCall<callPx;"OTM call < ATM call"];

/ Invalid inputs
badResult:@[{.commodity.black76.price[`call;x;75f;0.25;0.3;0.05]};-10f;{`ERROR}];
.testutil.assertTrue[badResult~`ERROR;"negative forward fails"];

-1 "PASS test_black76_price: call=",string[callPx],", put=",string putPx;
