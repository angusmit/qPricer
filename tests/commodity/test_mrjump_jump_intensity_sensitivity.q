\l core/init.q
/ Jump-intensity monotonicity. Increasing lambda holding everything else fixed
/ must:
/   (a) increase expectedJumpCount = lambda * T            (analytical, exact)
/   (b) increase spikeProbability  = 1 - exp(-lambda * T)  (analytical, exact)
/   (c) lift the MC call price when jumpMean is positive   (noisy; sized large)

x0:log 70f;
strikeVal:75f;
expiryVal:1f;
rateVal:0.05;
baseParams:`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility!(2f;log 75f;0.30;0.5;0.05;0.20);
highParams:@[baseParams;`jumpIntensity;:;3f];

baseExpectedJumps:.commodity.mrjump.expectedJumpCount[baseParams;expiryVal];
highExpectedJumps:.commodity.mrjump.expectedJumpCount[highParams;expiryVal];
.testutil.assertTrue[highExpectedJumps>baseExpectedJumps;"expectedJumpCount monotone in lambda"];

baseSpikeProb:.commodity.mrjump.spikeProbability[baseParams;expiryVal];
highSpikeProb:.commodity.mrjump.spikeProbability[highParams;expiryVal];
.testutil.assertTrue[highSpikeProb>baseSpikeProb;"spikeProbability monotone in lambda"];
.testutil.assertTrue[baseSpikeProb<1f;"baseSpikeProb in (0,1)"];
.testutil.assertTrue[highSpikeProb<1f;"highSpikeProb in (0,1)"];

mcConfig:.montecarlo.defaultMcConfig[];
mcConfig:@[mcConfig;`pathCount;:;30000];
mcConfig:@[mcConfig;`timeStepCount;:;50];

baseCallResult:.commodity.mrjump.europeanOptionPriceMC[`call;x0;strikeVal;expiryVal;rateVal;baseParams;mcConfig];
highCallResult:.commodity.mrjump.europeanOptionPriceMC[`call;x0;strikeVal;expiryVal;rateVal;highParams;mcConfig];
.testutil.assertTrue[highCallResult[`price]>baseCallResult`price;"MC call price monotone in lambda with positive jumpMean"];

-1 "PASS test_mrjump_jump_intensity_sensitivity: lowLambda=",string[baseParams`jumpIntensity],", highLambda=",string[highParams`jumpIntensity],", lowSpikeProb=",string[baseSpikeProb],", highSpikeProb=",string[highSpikeProb],", lowCall=",string[baseCallResult`price],", highCall=",string highCallResult`price;
