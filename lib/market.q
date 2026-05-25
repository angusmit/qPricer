/ market.q — flat market data: creation, access, bumping

.market.__requiredFields:`underlying`spot`riskFreeRate`dividendYield`volatility;

.market.createFlatMarketData:{[underlying;spot;riskFreeRate;dividendYield;volatility]
    md:`underlying`spot`riskFreeRate`dividendYield`volatility!(underlying;spot;riskFreeRate;dividendYield;volatility);
    .market.validateFlatMarketData md; md
 };

.market.validateFlatMarketData:{[md]
    .utilities.requireKeys[md;.market.__requiredFields;"marketData"];
    .utilities.assertPositive[md`spot;"spot"];
    .utilities.assertPositive[md`volatility;"volatility"];
    .utilities.assertNonNegative[md`dividendYield;"dividendYield"];
 };

/ Accessors — extra params ignored in v1, present for future API compatibility
.market.getSpot:{[md;underlying] md`spot};
.market.getRiskFreeRate:{[md;expiry] md`riskFreeRate};
.market.getDividendYield:{[md;underlying;expiry] md`dividendYield};
.market.getVolatility:{[md;underlying;strike;expiry] md`volatility};

/ Bumpers
.market.bumpSpot:{[md;underlying;bump] @[md;`spot;*;1f+bump]};
.market.bumpVolatility:{[md;bump] @[md;`volatility;+;bump]};
.market.bumpRiskFreeRate:{[md;bump] @[md;`riskFreeRate;+;bump]};
