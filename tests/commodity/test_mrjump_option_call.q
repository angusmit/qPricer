\l core/init.q
/ MC European call price is positive and has a meaningful confidence interval.
/ Higher jumpIntensity with positive jumpMean should lift the call price (positive
/ contribution to terminal price expectation via the jump compensator-free drift).

x0:log 70f;
strikeVal:75f;
expiryVal:1f;
rateVal:0.05;
baseParams:`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility!(2f;log 75f;0.30;1f;0.05;0.25);

mcConfig:.montecarlo.defaultMcConfig[];
mcConfig:@[mcConfig;`pathCount;:;20000];
mcConfig:@[mcConfig;`timeStepCount;:;50];

mcResult:.commodity.mrjump.europeanOptionPriceMC[`call;x0;strikeVal;expiryVal;rateVal;baseParams;mcConfig];
.testutil.assertTrue[mcResult[`price]>0f;"call price positive"];
.testutil.assertTrue[mcResult[`standardError]>0f;"standardError positive"];
.testutil.assertTrue[mcResult[`upperConfidence]>mcResult`lowerConfidence;"confidence band ordered"];

highIntensityParams:@[baseParams;`jumpIntensity;:;3f];
highResult:.commodity.mrjump.europeanOptionPriceMC[`call;x0;strikeVal;expiryVal;rateVal;highIntensityParams;mcConfig];
.testutil.assertTrue[highResult[`price]>mcResult`price;"higher intensity with positive jumpMean lifts call price"];

-1 "PASS test_mrjump_option_call: base ",string[mcResult`price]," (SE ",string[mcResult`standardError],"), high-intensity ",string highResult`price;
