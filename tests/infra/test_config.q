/ test_config.q - default configs and merge
\l lib/init.q

/ 1. Default pricing config has required fields
defaultCfg:.config.defaultPricingConfig[];
.testutil.assertEqual[defaultCfg`method;`explicit;"default method"];
.testutil.assertTrue[defaultCfg[`numberOfSpotSteps]>0;"default spot steps positive"];
.testutil.assertTrue[defaultCfg[`numberOfTimeSteps]>0;"default time steps positive"];

/ 2. Default passes validation
.config.validatePricingConfig defaultCfg;

/ 3. Default implied vol config
ivCfg:.config.defaultImpliedVolConfig[];
.testutil.assertTrue[ivCfg[`tolerance]>0f;"IV tolerance positive"];
.testutil.assertTrue[ivCfg[`maximumIterations]>0;"IV maxIter positive"];

/ 4. mergeConfig
merged:.config.mergeConfig[defaultCfg;`method`numberOfSpotSteps!(`crankNicolson;400)];
.testutil.assertEqual[merged`method;`crankNicolson;"merged method"];
.testutil.assertEqual[merged`numberOfSpotSteps;400;"merged spot steps"];
.testutil.assertEqual[merged`numberOfTimeSteps;defaultCfg`numberOfTimeSteps;"merged preserves unoverridden"];

/ 5. Merged config passes validation
.config.validatePricingConfig merged;

-1 "PASS test_config";
