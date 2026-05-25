/ product.q — option trade definition and validation

.product.__requiredFields:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional;
.product.__supportedProductTypes:enlist `equityOption;
.product.__supportedExerciseStyles:enlist `european;
.product.__supportedOptionTypes:`call`put;

.product.createOptionTrade:{[d]
    .product.validateOptionTrade d; d
 };

.product.validateOptionTrade:{[trade]
    .utilities.requireKeys[trade;.product.__requiredFields;"trade"];
    .utilities.assertSupportedValue[trade`productType;.product.__supportedProductTypes;"productType"];
    .utilities.assertSupportedValue[trade`exerciseStyle;.product.__supportedExerciseStyles;"exerciseStyle"];
    .utilities.assertSupportedValue[trade`optionType;.product.__supportedOptionTypes;"optionType"];
    .utilities.assertPositive[trade`strike;"strike"];
    .utilities.assertPositive[trade`expiry;"expiry"];
    .utilities.assertPositive[trade`notional;"notional"];
 };

.product.isCallOption:{[trade] trade[`optionType]~`call};
.product.isPutOption:{[trade] trade[`optionType]~`put};
.product.isEuropeanOption:{[trade] trade[`exerciseStyle]~`european};
