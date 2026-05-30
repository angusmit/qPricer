\l core/init.q
trade:.modelcheck.__defaultTrade[];
mkt:.modelcheck.__defaultMarket[];
mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(10000;50;42;0b;0b;0.95);
cfg:enlist[`mcConfig]!enlist mcConfig;

r1:.modelcheck.checkBatesHestonLimit[trade;mkt;cfg];
.testutil.assertTrue[r1[`status]~`OK;"Bates->Heston status OK"];
.testutil.assertTrue[`absoluteDifference in key r1;"has absoluteDifference"];
.testutil.assertTrue[`relativeDifference in key r1;"has relativeDifference"];

r2:.modelcheck.checkBatesMertonLimit[trade;mkt;cfg];
.testutil.assertTrue[r2[`status]~`OK;"Bates->Merton status OK"];

-1 "PASS test_model_limit_heston_merton_bates: batesHestonDiff=",string[r1`absoluteDifference],", batesMertonDiff=",string r2`absoluteDifference;
