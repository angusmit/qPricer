/ model.q — model definition

.model.createBlackScholesModel:{[]
    `modelName`volatilityType`rateType`dividendType!(`blackScholes;`flatVolatility;`flatRate;`flatDividend)
 };

.model.validateModel:{[model]
    .utilities.requireKeys[model;`modelName`volatilityType`rateType`dividendType;"model"];
    .utilities.assertSupportedValue[model`modelName;enlist `blackScholes;"modelName"];
    .utilities.assertSupportedValue[model`volatilityType;enlist `flatVolatility;"volatilityType"];
    .utilities.assertSupportedValue[model`rateType;enlist `flatRate;"rateType"];
    .utilities.assertSupportedValue[model`dividendType;enlist `flatDividend;"dividendType"];
 };
