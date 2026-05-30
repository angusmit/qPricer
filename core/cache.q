/ cache.q - pricing result cache (dict-based, explicit pass-through)

.cache.emptyCache:{[] ()!() };

.cache.createPricingKey:{[tradeRow;marketRow;fdmConfig]
    / Trade fields
    tradeKey:"|" sv string (tradeRow`underlying;tradeRow`optionType;tradeRow`exerciseStyle;
        tradeRow`strike;tradeRow`expiry;.product.getBarrierType tradeRow;.product.getBarrierLevel tradeRow);
    / Market fields (skip if no volatility — local-vol is not cacheable)
    if[not `volatility in key marketRow; :`$"NOCACHE"];
    marketKey:"|" sv string (marketRow`spot;marketRow`riskFreeRate;marketRow`dividendYield;marketRow`volatility);
    / Config fields
    configKey:"|" sv string (fdmConfig`method;fdmConfig`numberOfSpotSteps;fdmConfig`numberOfTimeSteps;
        fdmConfig`minimumSpot;fdmConfig`maximumSpot);
    `$tradeKey,"|",marketKey,"|",configKey
 };

.cache.get:{[cacheDict;pricingKey]
    cacheDict pricingKey
 };

.cache.put:{[cacheDict;pricingKey;pricingResult]
    cacheDict,enlist[pricingKey]!enlist pricingResult
 };

.cache.contains:{[cacheDict;pricingKey]
    pricingKey in key cacheDict
 };

.cache.getOrPrice:{[cacheDict;tradeRow;marketRow;pricingModel;fdmConfig]
    pricingKey:.cache.createPricingKey[tradeRow;marketRow;fdmConfig];
    / Not cacheable (local-vol) — price and return
    if[pricingKey~`NOCACHE;
        pricingResult:.engine.priceOption[tradeRow;marketRow;pricingModel;fdmConfig];
        :`cacheDict`pricingResult`cacheHit!(cacheDict;pricingResult;0b)];
    / Check cache
    if[.cache.contains[cacheDict;pricingKey];
        :`cacheDict`pricingResult`cacheHit!(cacheDict;.cache.get[cacheDict;pricingKey];1b)];
    / Price and store
    pricingResult:.engine.priceOption[tradeRow;marketRow;pricingModel;fdmConfig];
    updatedCache:.cache.put[cacheDict;pricingKey;pricingResult];
    `cacheDict`pricingResult`cacheHit!(updatedCache;pricingResult;0b)
 };

.cache.cacheSize:{[cacheDict] count cacheDict };
