/ model.q - model definition (Black-Scholes + local volatility)

.model.createBlackScholesModel:{[]
    `modelName`volatilityType`rateType`dividendType!(`blackScholes;`flatVolatility;`flatRate;`flatDividend)
 };

.model.createLocalVolatilityModel:{[]
    `modelName`volatilityType`rateType`dividendType!(`localVolatility;`localVolatility;`flatRate;`flatDividend)
 };

.model.validateModel:{[model]
    .utilities.requireKeys[model;`modelName`volatilityType`rateType`dividendType;"model"];
    .utilities.assertSupportedValue[model`modelName;`blackScholes`localVolatility;"modelName"];
    if[model[`modelName]~`blackScholes;
        .utilities.assertSupportedValue[model`volatilityType;enlist `flatVolatility;"volatilityType"]];
    if[model[`modelName]~`localVolatility;
        .utilities.assertSupportedValue[model`volatilityType;enlist `localVolatility;"volatilityType"]];
    .utilities.assertSupportedValue[model`rateType;enlist `flatRate;"rateType"];
    .utilities.assertSupportedValue[model`dividendType;enlist `flatDividend;"dividendType"];
 };

.model.isBlackScholes:{[model] model[`modelName]~`blackScholes};
.model.isLocalVolatility:{[model] model[`modelName]~`localVolatility};
