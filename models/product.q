/ product.q - option trade definition and validation (v0.14)

.product.__requiredFields:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional;
.product.__supportedProductTypes:enlist `equityOption;
.product.__supportedExerciseStyles:`european`american;
.product.__supportedOptionTypes:`call`put;
.product.__supportedBarrierTypes:`none`upAndOut`downAndOut`upAndIn`downAndIn;
.product.__knockOutTypes:`upAndOut`downAndOut;
.product.__knockInTypes:`upAndIn`downAndIn;

.product.createOptionTrade:{[tradeDictionary]
    .product.validateOptionTrade tradeDictionary; tradeDictionary
 };

.product.validateOptionTrade:{[trade]
    .utilities.requireKeys[trade;.product.__requiredFields;"trade"];
    .utilities.assertSupportedValue[trade`productType;.product.__supportedProductTypes;"productType"];
    .utilities.assertSupportedValue[trade`exerciseStyle;.product.__supportedExerciseStyles;"exerciseStyle"];
    .utilities.assertSupportedValue[trade`optionType;.product.__supportedOptionTypes;"optionType"];
    .utilities.assertPositive[trade`strike;"strike"];
    .utilities.assertPositive[trade`expiry;"expiry"];
    if[trade[`notional]=0f; '"notional must not be zero"];
    if[.product.isBarrierOption trade; .product.__validateBarrierFields trade];
 };

.product.__validateBarrierFields:{[trade]
    barrierType:.product.getBarrierType trade;
    .utilities.assertSupportedValue[barrierType;.product.__supportedBarrierTypes;"barrierType"];
    / Barrier options must be European
    if[not trade[`exerciseStyle]~`european;
        '"Barrier options require exerciseStyle = `european"];
    if[not `barrierLevel in key trade; '"Missing required field for barrier option: barrierLevel"];
    barrierLevel:trade`barrierLevel;
    .utilities.assertPositive[barrierLevel;"barrierLevel"];
    barrierRebate:.product.getBarrierRebate trade;
    if[barrierRebate<>0f; '"Non-zero barrier rebate is not supported"];
 };

.product.isCallOption:{[trade] trade[`optionType]~`call};
.product.isPutOption:{[trade] trade[`optionType]~`put};
.product.isEuropeanOption:{[trade] trade[`exerciseStyle]~`european};
.product.isAmericanOption:{[trade] trade[`exerciseStyle]~`american};

.product.isBarrierOption:{[trade]
    barrierType:.product.getBarrierType trade;
    not barrierType~`none
 };

.product.isKnockIn:{[trade]
    (.product.getBarrierType trade) in .product.__knockInTypes
 };

.product.isKnockOut:{[trade]
    (.product.getBarrierType trade) in .product.__knockOutTypes
 };

.product.knockInToKnockOut:{[barrierType]
    $[barrierType~`upAndIn;`upAndOut;
      barrierType~`downAndIn;`downAndOut;
      '"Not a knock-in type: ",string barrierType]
 };

.product.getBarrierType:{[trade]
    if[not `barrierType in key trade; :`none];
    trade`barrierType
 };

.product.getBarrierLevel:{[trade]
    if[not `barrierLevel in key trade; :0Nf];
    trade`barrierLevel
 };

.product.getBarrierRebate:{[trade]
    if[not `rebate in key trade; :0f];
    trade`rebate
 };
