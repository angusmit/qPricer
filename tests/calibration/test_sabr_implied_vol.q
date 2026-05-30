/ test_sabr_implied_vol.q
\l core/init.q
sabrParams:`alpha`beta`rho`nu!(0.2;0.5;-0.3;0.4);

/ ATM
atmVol:.sabr.impliedVolAtm[100f;1f;sabrParams];
.testutil.assertTrue[atmVol>0f;"ATM vol positive"];
.testutil.assertTrue[not null atmVol;"ATM vol finite"];

/ Non-ATM
otmVol:.sabr.impliedVol[100f;120f;1f;sabrParams];
.testutil.assertTrue[otmVol>0f;"OTM vol positive"];
itmVol:.sabr.impliedVol[100f;80f;1f;sabrParams];
.testutil.assertTrue[itmVol>0f;"ITM vol positive"];

/ Near-ATM stability
nearAtmVol:.sabr.impliedVol[100f;100.0001;1f;sabrParams];
.testutil.assertTrue[nearAtmVol>0f;"near-ATM vol positive"];
.testutil.assertNear[nearAtmVol;atmVol;0.001;"near-ATM close to ATM"];

/ Increasing nu changes smile
lowNuParams:`alpha`beta`rho`nu!(0.2;0.5;-0.3;0.1);
highNuParams:`alpha`beta`rho`nu!(0.2;0.5;-0.3;0.8);
lowNuVol:.sabr.impliedVol[100f;120f;1f;lowNuParams];
highNuVol:.sabr.impliedVol[100f;120f;1f;highNuParams];
.testutil.assertTrue[not lowNuVol=highNuVol;"nu changes smile"];

-1 "PASS test_sabr_implied_vol: ATM=",string[atmVol],", OTM=",string[otmVol],", ITM=",string itmVol;
