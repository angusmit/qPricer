\l lib/init.q
trade:.modelcheck.__defaultTrade[];
mkt:.modelcheck.__defaultMarket[];

r1:.modelcheck.checkLocalVolFlatLimit[trade;mkt;()!()];
.testutil.assertTrue[r1`passed;"Local vol flat limit passes"];
.testutil.assertTrue[r1[`absoluteDifference]<0.001;"local vol diff tiny"];

r2:.modelcheck.checkSabrFlatSmileLimit[()!();()!()];
.testutil.assertTrue[r2[`status]~`OK;"SABR flat status OK"];
.testutil.assertTrue[r2[`absoluteDifference]<0.01;"SABR near-flat smile variation small"];

-1 "PASS test_model_limit_local_vol_sabr: lvDiff=",string[r1`absoluteDifference],", sabrDiff=",string r2`absoluteDifference;
