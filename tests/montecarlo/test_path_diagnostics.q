/ test_path_diagnostics.q
\l core/init.q
mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(1000;50;42;0b;0b;0.95);
pathMatrix:.montecarlo.simulateGBMPaths[100f;0.05;0f;0.2;1f;mcConfig];

minVec:.pathdiag.pathMinimum pathMatrix;
maxVec:.pathdiag.pathMaximum pathMatrix;
finalVec:.pathdiag.pathFinal pathMatrix;
avgVec:.pathdiag.pathAverage pathMatrix;

.testutil.assertTrue[1000=count minVec;"1000 min values"];
.testutil.assertTrue[1000=count maxVec;"1000 max values"];
.testutil.assertTrue[1000=count finalVec;"1000 final values"];
.testutil.assertTrue[all maxVec>=minVec;"max >= min for all paths"];
.testutil.assertTrue[all finalVec>0f;"all final spots positive"];

summaryResult:.pathdiag.summary[pathMatrix;1f];
.testutil.assertTrue[summaryResult[`pathCount]=1000;"summary pathCount"];
.testutil.assertTrue[summaryResult[`timeStepCount]=50;"summary timeStepCount"];
.testutil.assertTrue[summaryResult[`averageFinal]>0f;"summary avgFinal positive"];

-1 "PASS test_path_diagnostics: avgFinal=",string[summaryResult`averageFinal],", avgMin=",string[summaryResult`averageMinimum],", avgMax=",string summaryResult`averageMaximum;
