\l lib/init.q
/ Closed-form futures price F(0,T) should match Monte Carlo mean of exp(X_T).

params:`meanReversionSpeed`longRunLogMean`volatility!(1.5;log 75f;0.30);
x0:log 70f;
tauVal:1f;

mcConfig:.montecarlo.defaultMcConfig[];
mcConfig:@[mcConfig;`pathCount;:;30000];
mcConfig:@[mcConfig;`timeStepCount;:;1];

pathMatrix:.commodity.schwartz.simulatePaths[x0;params;tauVal;mcConfig];
terminalSpots:last each pathMatrix;
mcMean:avg terminalSpots;
mcStdErr:(dev terminalSpots)%sqrt count terminalSpots;
analyticalFwd:.commodity.schwartz.futuresPrice[x0;params;tauVal];

.testutil.assertTrue[analyticalFwd>0f;"forward positive"];
.testutil.assertTrue[(abs analyticalFwd-mcMean)<4f*mcStdErr;"closed-form forward within 4 SE of MC mean"];

-1 "PASS test_schwartz_futures_price: closed ",string[analyticalFwd]," vs MC ",string[mcMean]," (SE ",string[mcStdErr],")";
