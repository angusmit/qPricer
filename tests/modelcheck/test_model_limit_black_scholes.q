\l lib/init.q
trade:.modelcheck.__defaultTrade[];
mkt:.modelcheck.__defaultMarket[];
mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(10000;50;42;0b;0b;0.95);
cfg:enlist[`mcConfig]!enlist mcConfig;

r1:.modelcheck.checkHestonBlackScholesLimit[trade;mkt;cfg];
.testutil.assertTrue[r1`passed;"Heston->BS limit passes"];
.testutil.assertTrue[r1[`status]~`OK;"Heston->BS status OK"];

r2:.modelcheck.checkMertonBlackScholesLimit[trade;mkt;cfg];
.testutil.assertTrue[r2`passed;"Merton->BS limit passes"];

r3:.modelcheck.checkBatesBlackScholesLimit[trade;mkt;cfg];
.testutil.assertTrue[r3`passed;"Bates->BS limit passes"];

/ Check required columns
firstRow:r1;
.testutil.assertTrue[`checkName in key firstRow;"has checkName"];
.testutil.assertTrue[`passed in key firstRow;"has passed"];
.testutil.assertTrue[`status in key firstRow;"has status"];

-1 "PASS test_model_limit_black_scholes: hestonDiff=",string[r1`absoluteDifference],", mertonDiff=",string[r2`absoluteDifference],", batesDiff=",string r3`absoluteDifference;
