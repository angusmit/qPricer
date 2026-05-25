/ config.q - finite-difference configuration

.config.__requiredFields:`method`numberOfSpotSteps`numberOfTimeSteps`maximumSpot;
.config.__supportedMethods:`explicit`crankNicolson;
.config.__supportedInterpolation:enlist `linear;

.config.createFiniteDifferenceConfig:{[configDict]
    cfg:.config.withDefaults configDict;
    .config.validateFiniteDifferenceConfig cfg; cfg
 };

.config.withDefaults:{[cfg]
    defaults:`minimumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(0f;`linear;1b;1b);
    missingKeys:(key defaults) where not (key defaults) in key cfg;
    cfg,missingKeys#defaults
 };

.config.validateFiniteDifferenceConfig:{[cfg]
    .utilities.requireKeys[cfg;.config.__requiredFields;"config"];
    .utilities.assertSupportedValue[cfg`method;.config.__supportedMethods;"method"];
    .utilities.assertPositiveInteger[cfg`numberOfSpotSteps;"numberOfSpotSteps"];
    .utilities.assertPositiveInteger[cfg`numberOfTimeSteps;"numberOfTimeSteps"];
    .utilities.assertNonNegative[cfg`minimumSpot;"minimumSpot"];
    .utilities.assertPositive[cfg`maximumSpot;"maximumSpot"];
    if[cfg[`maximumSpot]<=cfg`minimumSpot; '"maximumSpot must be greater than minimumSpot"];
    .utilities.assertSupportedValue[cfg`interpolationMethod;.config.__supportedInterpolation;"interpolationMethod"];
    .utilities.assertBoolean[cfg`returnFullGrid;"returnFullGrid"];
    .utilities.assertBoolean[cfg`stabilityCheck;"stabilityCheck"];
 };
