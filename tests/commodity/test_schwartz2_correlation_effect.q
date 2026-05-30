\l core/init.q
/ Correlation increases log-price variance via the cross-covariance term, hence
/ increases European option price (vega positive). Sweep rho and verify strict
/ ordering. Also verify validateParams rejects rho at or outside +/-1.

shortFactor0:0f;
longFactor0:log 75f;
strikeVal:75f;
expiryVal:1f;
rateVal:0.05;
baseParams:`meanReversionSpeed`shortVolatility`longVolatility`longDrift`correlation!(1.5;0.35;0.15;0f;0f);

paramsNeg:@[baseParams;`correlation;:;-0.5];
paramsZero:baseParams;
paramsPos:@[baseParams;`correlation;:;0.5];

callNeg:.commodity.schwartz2.europeanOptionPrice[`call;shortFactor0;longFactor0;strikeVal;expiryVal;rateVal;paramsNeg];
callZero:.commodity.schwartz2.europeanOptionPrice[`call;shortFactor0;longFactor0;strikeVal;expiryVal;rateVal;paramsZero];
callPos:.commodity.schwartz2.europeanOptionPrice[`call;shortFactor0;longFactor0;strikeVal;expiryVal;rateVal;paramsPos];

.testutil.assertTrue[callNeg<callZero;"call increases when correlation rises from -0.5 to 0"];
.testutil.assertTrue[callZero<callPos;"call increases when correlation rises from 0 to 0.5"];

paramsRhoLow:@[baseParams;`correlation;:;-1f];
rLow:@[.commodity.schwartz2.validateParams;paramsRhoLow;{`ERROR}];
.testutil.assertTrue[rLow~`ERROR;"correlation = -1 rejected"];

paramsRhoHigh:@[baseParams;`correlation;:;1f];
rHigh:@[.commodity.schwartz2.validateParams;paramsRhoHigh;{`ERROR}];
.testutil.assertTrue[rHigh~`ERROR;"correlation = +1 rejected"];

paramsRhoOver:@[baseParams;`correlation;:;1.5];
rOver:@[.commodity.schwartz2.validateParams;paramsRhoOver;{`ERROR}];
.testutil.assertTrue[rOver~`ERROR;"correlation > 1 rejected"];

-1 "PASS test_schwartz2_correlation_effect: calls -0.5,0,+0.5 = ",string[callNeg]," ",string[callZero]," ",string callPos;
