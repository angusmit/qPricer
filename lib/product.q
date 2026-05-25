/ product.q — option trade definition and validation

.product.__requiredFields:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional;
.product.__supportedProductTypes:enlist `equityOption;
.product.__supportedExerciseStyles:`european`american;
.product.__supportedOptionTypes:`call`put;

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
    / American call on non-dividend stock = European call; allow but warn
    if[(trade[`exerciseStyle]~`american) and trade[`optionType]~`call;
        -1 "WARNING: American call on non-dividend stock has no early exercise premium"];
 };

.product.isCallOption:{[trade] trade[`optionType]~`call};
.product.isPutOption:{[trade] trade[`optionType]~`put};
.product.isEuropeanOption:{[trade] trade[`exerciseStyle]~`european};
.product.isAmericanOption:{[trade] trade[`exerciseStyle]~`american};
