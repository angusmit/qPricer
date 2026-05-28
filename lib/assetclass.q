/ assetclass.q - asset class registry and routing (v0.35)
/ Registered commodity model families:
/   black76            - Black-76 futures option
/   kirk               - Kirk spread approximation
/   meanReversion      - generic mean-reverting marker
/   schwartzOneFactor  - .commodity.schwartz       (v0.33, log mean reversion)
/   schwartzTwoFactor  - .commodity.schwartz2      (v0.34, spot + stochastic convenience yield)
/   meanRevertingJump  - .commodity.mrjump         (v0.35, OU log-price + compound-Poisson jumps)

.assetclass.__productMap:`equityOption`asianOption`basketOption`lookbackOption`commodityFuture`commodityOption`spreadOption`sparkSpreadOption`powerForward!`equity`equity`equity`equity`commodity`commodity`commodity`electricity`electricity;
.assetclass.__modelMap:`blackScholes`localVolatility`heston`merton`bates`sabr`black76`kirk`meanReversion`schwartzOneFactor`schwartzTwoFactor`meanRevertingJump!`equity`equity`equity`equity`equity`equity`commodity`commodity`commodity`commodity`commodity`commodity;
.assetclass.__supported:`equity`commodity`electricity`rates`fx;

.assetclass.supportedAssetClasses:{[] .assetclass.__supported};

.assetclass.assetClassForProduct:{[productType]
    if[productType in key .assetclass.__productMap; :.assetclass.__productMap productType];
    `unknown
 };

.assetclass.modelFamily:{[modelType]
    if[modelType in key .assetclass.__modelMap; :.assetclass.__modelMap modelType];
    `unknown
 };

.assetclass.validateAssetClass:{[assetClass]
    if[not assetClass in .assetclass.__supported; '"Unsupported asset class: ",string assetClass];
 };

.assetclass.routingHint:{[productType;modelType]
    `assetClass`productFamily`modelFamily!(.assetclass.assetClassForProduct productType;productType;.assetclass.modelFamily modelType)
 };
