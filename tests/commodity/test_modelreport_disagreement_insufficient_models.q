\l lib/init.q
/ Status policy:
/   zero OK rows                 -> `ERROR
/   1 .. minimumOkModels-1 rows  -> `warning
/   minimumOkModels or more      -> `OK
/ Exercise all three branches by synthesising controlled priceAllModels-shaped
/ tables.

priceRowOk:`modelName`modelPrice`standardError`lowerConfidence`upperConfidence`pricingStatus`pricingErrorMessage!(
    `black76;5.0;0Nf;0Nf;0Nf;`OK;"");
priceRowFail:`modelName`modelPrice`standardError`lowerConfidence`upperConfidence`pricingStatus`pricingErrorMessage!(
    `schwartz;0Nf;0Nf;0Nf;0Nf;`ERROR;"missing");

singleOkTbl:(enlist priceRowOk),enlist priceRowFail;
allFailTbl:(enlist priceRowFail),enlist @[priceRowFail;`modelName;:;`schwartz2];

disCfg:.commodity.modelreport.defaultDisagreementConfig[];

okPriceTbl:(enlist priceRowOk),enlist @[priceRowOk;(`modelName;`modelPrice);:;(`schwartz;6.0)];
fullDis:.commodity.modelreport.priceDisagreement[okPriceTbl;disCfg];
.testutil.assertTrue[`OK=fullDis`status;"two OK rows at minimumOkModels=2 -> OK"];

partialDis:.commodity.modelreport.priceDisagreement[singleOkTbl;disCfg];
.testutil.assertTrue[`warning=partialDis`status;"one OK row below minimumOkModels -> warning"];
.testutil.assertTrue[1=partialDis`okModelCount;"partialDis okModelCount = 1"];

emptyDis:.commodity.modelreport.priceDisagreement[allFailTbl;disCfg];
.testutil.assertTrue[`ERROR=emptyDis`status;"zero OK rows -> ERROR"];
.testutil.assertTrue[0=emptyDis`okModelCount;"emptyDis okModelCount = 0"];

scenRowOk:`modelName`basePrice`scenarioPrice`scenarioPnL`quantity`contractMultiplier`status`errorMessage!(
    `black76;5.0;5.5;500f;1f;1000f;`OK;"");
scenRowFail:`modelName`basePrice`scenarioPrice`scenarioPnL`quantity`contractMultiplier`status`errorMessage!(
    `schwartz;0Nf;0Nf;0Nf;1f;1000f;`ERROR;"");
singleOkScen:(enlist scenRowOk),enlist scenRowFail;
allFailScen:(enlist scenRowFail),enlist @[scenRowFail;`modelName;:;`schwartz2];

partialScen:.commodity.modelreport.scenarioDisagreement[singleOkScen;disCfg];
.testutil.assertTrue[`warning=partialScen`status;"scenarioDisagreement partial -> warning"];
emptyScen:.commodity.modelreport.scenarioDisagreement[allFailScen;disCfg];
.testutil.assertTrue[`ERROR=emptyScen`status;"scenarioDisagreement empty -> ERROR"];

-1 "PASS test_modelreport_disagreement_insufficient_models: warningStatus=",string[partialDis`status],", errorStatus=",string[emptyDis`status];
