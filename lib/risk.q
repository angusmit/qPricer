/ risk.q — scenario risk reporting
/ Bumps market data and reprices through the public engine API.
/ Conventions:
/   spot bumps: relative (0.01 = +1%)
/   vol bumps:  absolute (0.01 = +1 vol point)
/   rate bumps: absolute (0.0025 = +25bp)

/ --- Public ---

.risk.generateScenarioReport:{[trade;marketData;model;config]
    .product.validateOptionTrade trade;
    .market.validateFlatMarketData marketData;
    .model.validateModel model;
    .config.validateFiniteDifferenceConfig config;
    basePriceResult:.engine.priceOption[trade;marketData;model;config];
    baseUnitPrice:basePriceResult`unitPrice;
    baseNotionalPrice:basePriceResult`notionalPrice;
    scenarioDefinitions:.risk.__createScenarioDefinitions[];
    scenarioRows:.risk.__priceScenario[trade;marketData;model;config;baseUnitPrice;baseNotionalPrice;] each scenarioDefinitions;
    ([] tradeId:scenarioRows[;`tradeId];
        underlying:scenarioRows[;`underlying];
        optionType:scenarioRows[;`optionType];
        scenario:scenarioRows[;`scenario];
        spot:scenarioRows[;`spot];
        riskFreeRate:scenarioRows[;`riskFreeRate];
        dividendYield:scenarioRows[;`dividendYield];
        volatility:scenarioRows[;`volatility];
        unitPrice:scenarioRows[;`unitPrice];
        unitPnL:scenarioRows[;`unitPnL];
        notionalPrice:scenarioRows[;`notionalPrice];
        notionalPnL:scenarioRows[;`notionalPnL])
 };

/ --- Internal ---

.risk.__createScenarioDefinitions:{[]
    (
        `scenario`spotBump`volatilityBump`rateBump!(`base;0f;0f;0f);
        `scenario`spotBump`volatilityBump`rateBump!(`spotUp1Pct;0.01;0f;0f);
        `scenario`spotBump`volatilityBump`rateBump!(`spotDown1Pct;-0.01;0f;0f);
        `scenario`spotBump`volatilityBump`rateBump!(`spotUp5Pct;0.05;0f;0f);
        `scenario`spotBump`volatilityBump`rateBump!(`spotDown5Pct;-0.05;0f;0f);
        `scenario`spotBump`volatilityBump`rateBump!(`volatilityUp1Point;0f;0.01;0f);
        `scenario`spotBump`volatilityBump`rateBump!(`volatilityDown1Point;0f;-0.01;0f);
        `scenario`spotBump`volatilityBump`rateBump!(`rateUp25Bp;0f;0f;0.0025);
        `scenario`spotBump`volatilityBump`rateBump!(`rateDown25Bp;0f;0f;-0.0025)
    )
 };

.risk.__applyScenario:{[marketData;scenarioDefinition]
    scenarioMkt:marketData;
    if[0<>scenarioDefinition`spotBump;
        scenarioMkt:.market.bumpSpot[scenarioMkt;scenarioMkt`underlying;scenarioDefinition`spotBump]];
    if[0<>scenarioDefinition`volatilityBump;
        scenarioMkt:.market.bumpVolatility[scenarioMkt;scenarioDefinition`volatilityBump]];
    if[0<>scenarioDefinition`rateBump;
        scenarioMkt:.market.bumpRiskFreeRate[scenarioMkt;scenarioDefinition`rateBump]];
    scenarioMkt
 };

.risk.__priceScenario:{[trade;marketData;model;config;baseUnitPrice;baseNotionalPrice;scenarioDefinition]
    scenarioMkt:.risk.__applyScenario[marketData;scenarioDefinition];
    scenarioPriceResult:.engine.priceOption[trade;scenarioMkt;model;config];
    scenarioUnitPrice:scenarioPriceResult`unitPrice;
    scenarioNotionalPrice:scenarioPriceResult`notionalPrice;
    unitPnL:scenarioUnitPrice - baseUnitPrice;
    notionalPnL:scenarioNotionalPrice - baseNotionalPrice;
    `tradeId`underlying`optionType`scenario`spot`riskFreeRate`dividendYield`volatility`unitPrice`unitPnL`notionalPrice`notionalPnL!(
        trade`tradeId;
        trade`underlying;
        trade`optionType;
        scenarioDefinition`scenario;
        scenarioMkt`spot;
        scenarioMkt`riskFreeRate;
        scenarioMkt`dividendYield;
        scenarioMkt`volatility;
        scenarioUnitPrice;
        unitPnL;
        scenarioNotionalPrice;
        notionalPnL)
 };
