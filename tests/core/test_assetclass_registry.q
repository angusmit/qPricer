\l lib/init.q
supported:.assetclass.supportedAssetClasses[];
.testutil.assertTrue[`equity in supported;"has equity"];
.testutil.assertTrue[`commodity in supported;"has commodity"];
.testutil.assertTrue[`electricity in supported;"has electricity"];

/ Product mapping
.testutil.assertTrue[`equity=.assetclass.assetClassForProduct `equityOption;"equityOption -> equity"];
.testutil.assertTrue[`commodity=.assetclass.assetClassForProduct `commodityFuture;"commodityFuture -> commodity"];
.testutil.assertTrue[`commodity=.assetclass.assetClassForProduct `commodityOption;"commodityOption -> commodity"];
.testutil.assertTrue[`commodity=.assetclass.assetClassForProduct `spreadOption;"spreadOption -> commodity"];

/ Model family
.testutil.assertTrue[`commodity=.assetclass.modelFamily `black76;"black76 -> commodity"];
.testutil.assertTrue[`equity=.assetclass.modelFamily `blackScholes;"blackScholes -> equity"];

/ Invalid
badAc:@[.assetclass.validateAssetClass;`crypto;{`ERROR}];
.testutil.assertTrue[badAc~`ERROR;"unknown asset class fails"];

/ Routing
hint:.assetclass.routingHint[`commodityOption;`black76];
.testutil.assertTrue[hint[`assetClass]=`commodity;"routing commodity"];

-1 "PASS test_assetclass_registry";
