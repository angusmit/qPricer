/ product.q - option trade definition and validation
/ Supports vanilla + barrier options (v0.5)

.product.__requiredFields:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional;
.product.__supportedProductTypes:enlist `equityOption;
.product.__supportedExerciseStyles:`european`american;
.product.__supportedOptionTypes:`call`put;
.product.__supportedBarrierTypes:`none`upAndOut`downAndOut;

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
    .utilities.assertPositive[trade`notional;"notional"];
    / Validate barrier fields if present
    if[.product.isBarrierOption trade; .product.__validateBarrierFields trade];
 };

.product.__validateBarrierFields:{[trade]
    barrierType:.product.getBarrierType trade;
    .utilities.assertSupportedValue[barrierType;.product.__supportedBarrierTypes;"barrierType"];
    / Must be European for barrier v0.5
    if[not trade[`exerciseStyle]~`european;
        '"Barrier options require exerciseStyle = `european in v0.5"];
    / Validate barrierLevel
    if[not `barrierLevel in key trade; '"Missing required field for barrier option: barrierLevel"];
    barrierLevel:trade`barrierLevel;
    .utilities.assertPositive[barrierLevel;"barrierLevel"];
    / Validate rebate
    barrierRebate:.product.getBarrierRebate trade;
    if[barrierRebate<>0f; '"Non-zero barrier rebate is not supported in v0.5"];
    / Validate combination
    if[(barrierType~`upAndOut) and not trade[`optionType]~`call;
        '"upAndOut barrier is only supported with call in v0.5"];
    if[(barrierType~`downAndOut) and not trade[`optionType]~`put;
        '"downAndOut barrier is only supported with put in v0.5"];
 };

.product.isCallOption:{[trade] trade[`optionType]~`call};
.product.isPutOption:{[trade] trade[`optionType]~`put};
.product.isEuropeanOption:{[trade] trade[`exerciseStyle]~`european};
.product.isAmericanOption:{[trade] trade[`exerciseStyle]~`american};

.product.isBarrierOption:{[trade]
    barrierType:.product.getBarrierType trade;
    not barrierType~`none
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
