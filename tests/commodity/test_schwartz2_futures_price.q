\l core/init.q
/ Closed-form futures price F(0,T) should match MC mean of exp(X+Y) within MC SE.

params:`meanReversionSpeed`shortVolatility`longVolatility`longDrift`correlation!(1.5;0.35;0.15;0f;0.3);
shortFactor0:0f;
longFactor0:log 75f;
tauVal:1f;

mcConfig:.montecarlo.defaultMcConfig[];
mcConfig:@[mcConfig;`pathCount;:;30000];
mcConfig:@[mcConfig;`timeStepCount;:;1];

pathMatrix:.commodity.schwartz2.simulatePaths[shortFactor0;longFactor0;params;tauVal;mcConfig];
terminalSpots:last each pathMatrix;
mcMean:avg terminalSpots;
mcStdErr:(dev terminalSpots)%sqrt count terminalSpots;
analyticalFwd:.commodity.schwartz2.futuresPrice[shortFactor0;longFactor0;params;tauVal];

.testutil.assertTrue[analyticalFwd>0f;"forward positive"];
.testutil.assertTrue[(abs analyticalFwd-mcMean)<4f*mcStdErr;"closed-form forward within 4 SE of MC mean"];

curveValues:.commodity.schwartz2.futuresCurve[shortFactor0;longFactor0;params;0.25 0.5 1 2 5f];
.testutil.assertTrue[all curveValues>0f;"futures curve strictly positive"];

-1 "PASS test_schwartz2_futures_price: closed ",string[analyticalFwd]," vs MC ",string[mcMean]," (SE ",string[mcStdErr],")";
