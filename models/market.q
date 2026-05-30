/ market.q - market data: creation, access, bumping
/ Supports flat (BS) and local volatility market data

.market.__flatRequiredFields:`underlying`spot`riskFreeRate`dividendYield`volatility;
.market.__localVolRequiredFields:`underlying`spot`riskFreeRate`dividendYield`localVolatilityFunction;

/ --- Creation ---

.market.createFlatMarketData:{[underlying;spot;riskFreeRate;dividendYield;volatility]
    flatMkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(underlying;spot;riskFreeRate;dividendYield;volatility);
    .market.validateFlatMarketData flatMkt; flatMkt
 };

.market.createLocalVolatilityMarketData:{[underlying;spot;riskFreeRate;dividendYield;localVolFunction]
    localMkt:`underlying`spot`riskFreeRate`dividendYield`localVolatilityFunction!(
        underlying;spot;riskFreeRate;dividendYield;localVolFunction);
    .market.validateLocalVolatilityMarketData localMkt; localMkt
 };

/ --- Validation ---

.market.validateFlatMarketData:{[marketData]
    .utilities.requireKeys[marketData;.market.__flatRequiredFields;"marketData"];
    .utilities.assertPositive[marketData`spot;"spot"];
    .utilities.assertPositive[marketData`volatility;"volatility"];
    .utilities.assertNonNegative[marketData`dividendYield;"dividendYield"];
 };

.market.validateLocalVolatilityMarketData:{[marketData]
    .utilities.requireKeys[marketData;.market.__localVolRequiredFields;"localVolatilityMarketData"];
    .utilities.assertPositive[marketData`spot;"spot"];
    .utilities.assertNonNegative[marketData`dividendYield;"dividendYield"];
    / Test that the function is callable and returns positive vol
    testVol:marketData[`localVolatilityFunction][marketData`spot;0f];
    if[not testVol>0f; '"localVolatilityFunction must return positive volatility"];
 };

.market.validateMarketData:{[marketData;model]
    / Vol surface market data is valid for BS model
    if[(`marketDataType in key marketData) and marketData[`marketDataType]~`volSurface;
        :.market.validateSurfaceMarketData marketData];
    if[model[`modelName]~`blackScholes; :.market.validateFlatMarketData marketData];
    if[model[`modelName]~`localVolatility; :.market.validateLocalVolatilityMarketData marketData];
    '"Unsupported model for market data validation: ",string model`modelName
 };

.market.validateSurfaceMarketData:{[marketData]
    .utilities.requireKeys[marketData;`underlying`spot`riskFreeRate`dividendYield`volSurface;"surfaceMarketData"];
    .utilities.assertPositive[marketData`spot;"spot"];
    .utilities.assertNonNegative[marketData`dividendYield;"dividendYield"];
    .surface.validateVolSurface marketData`volSurface;
 };

/ --- Accessors (extra params for future API compat) ---

.market.getSpot:{[marketData;underlying] marketData`spot};
.market.getRiskFreeRate:{[marketData;expiry] marketData`riskFreeRate};
.market.getDividendYield:{[marketData;underlying;expiry] marketData`dividendYield};
.market.getVolatility:{[marketData;underlying;strike;expiry]
    if[(`marketDataType in key marketData) and marketData[`marketDataType]~`volSurface;
        :.surface.getSurfaceVolatility[marketData`volSurface;strike;expiry]];
    marketData`volatility
 };

/ Local vol accessor: returns volatility vector for spot vector at given time
.market.getLocalVolatility:{[marketData;spotValues;timePoint]
    if[not `localVolatilityFunction in key marketData;
        '"Local volatility function not found in marketData"];
    marketData[`localVolatilityFunction][spotValues;timePoint]
 };

/ --- Bumpers ---

.market.bumpSpot:{[marketData;underlying;bumpFraction] @[marketData;`spot;*;1f+bumpFraction]};
.market.bumpVolatility:{[marketData;absoluteBump] @[marketData;`volatility;+;absoluteBump]};
.market.bumpRiskFreeRate:{[marketData;absoluteBump] @[marketData;`riskFreeRate;+;absoluteBump]};
