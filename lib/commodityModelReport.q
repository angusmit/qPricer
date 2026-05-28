/ commodityModelReport.q - cross-model commodity option comparison report (v0.36)
/ Reporting layer over .commodity.black76, .commodity.schwartz, .commodity.schwartz2
/ and .commodity.mrjump. Prices one option setup with each model, summarises
/ differences against a chosen baseline, applies simple scenario shifts, and
/ computes per-model PnL. This module changes no pricing formulas.
/ Public:  validateOptionSetup, priceBlack76, priceSchwartz, priceSchwartz2,
/          priceMrJump, priceAllModels, comparisonSummary, priceDifferences,
/          scenarioShift, scenarioPnL, runComparison, modelRanking.
/ Private: __priceRow, __failRow, __dictsToTable, __bumpModelField,
/          __bumpModelSubField, __makeAddFn, __makeScaleFn, __makeMulFn.

.commodity.modelreport.__supportedOptionTypes:`call`put;

/ --- Helpers --------------------------------------------------------

/ Common-shape row dict used by every pricing wrapper.
.commodity.modelreport.__priceRow:{[modelName;modelPriceVal;standardErrorVal;lowerConfVal;upperConfVal;pricingStatus;pricingErrorMessage]
    `modelName`modelPrice`standardError`lowerConfidence`upperConfidence`pricingStatus`pricingErrorMessage!(
        modelName;modelPriceVal;standardErrorVal;lowerConfVal;upperConfVal;pricingStatus;pricingErrorMessage)
 };

.commodity.modelreport.__failRow:{[modelName;errorMessage]
    .commodity.modelreport.__priceRow[modelName;0Nf;0Nf;0Nf;0Nf;`ERROR;errorMessage]
 };

/ Build a table from a list of dicts with identical key order.
.commodity.modelreport.__dictsToTable:{[rowDicts]
    if[0=count rowDicts; :()];
    (,/) enlist each rowDicts
 };

/ Bump a top-level field on inputs[modelKey] if the key chain exists.
/ bumpFn is a monadic function applied to the current scalar value.
.commodity.modelreport.__bumpModelField:{[inputs;modelKey;fieldKey;bumpFn]
    if[not modelKey in key inputs; :inputs];
    inputDict:inputs modelKey;
    if[not fieldKey in key inputDict; :inputs];
    newVal:bumpFn inputDict fieldKey;
    newInputDict:@[inputDict;fieldKey;:;newVal];
    @[inputs;modelKey;:;newInputDict]
 };

/ Bump a nested field on inputs[modelKey][outerKey][innerKey] if the chain exists.
.commodity.modelreport.__bumpModelSubField:{[inputs;modelKey;outerKey;innerKey;bumpFn]
    if[not modelKey in key inputs; :inputs];
    inputDict:inputs modelKey;
    if[not outerKey in key inputDict; :inputs];
    outerDict:inputDict outerKey;
    if[not innerKey in key outerDict; :inputs];
    newVal:bumpFn outerDict innerKey;
    newOuterDict:@[outerDict;innerKey;:;newVal];
    newInputDict:@[inputDict;outerKey;:;newOuterDict];
    @[inputs;modelKey;:;newInputDict]
 };

/ Closures for scenario bumps.
.commodity.modelreport.__makeAddFn:{[bumpAmount] {[bump;currentVal] currentVal+bump}[bumpAmount;]};
.commodity.modelreport.__makeScaleFn:{[shiftPct] {[pct;currentVal] currentVal*1f+pct}[shiftPct;]};
.commodity.modelreport.__makeMulFn:{[factor] {[m;currentVal] currentVal*m}[factor;]};

/ --- Validation -----------------------------------------------------

.commodity.modelreport.validateOptionSetup:{[optionSetup]
    requiredKeys:`optionType`strikePrice`expiry`riskFreeRate`forwardPrice`spotPrice;
    presentKeys:key optionSetup;
    missingKeys:requiredKeys where not requiredKeys in presentKeys;
    if[0<count missingKeys;
        '"modelreport optionSetup missing: ",", " sv string missingKeys];
    if[not optionSetup[`optionType] in .commodity.modelreport.__supportedOptionTypes;
        '"modelreport optionType must be call or put"];
    if[optionSetup[`strikePrice]<=0f; '"modelreport strikePrice must be positive"];
    if[optionSetup[`expiry]<=0f; '"modelreport expiry must be positive"];
    if[optionSetup[`forwardPrice]<=0f; '"modelreport forwardPrice must be positive"];
    if[optionSetup[`spotPrice]<=0f; '"modelreport spotPrice must be positive"];
 };

/ --- Per-model pricing wrappers ------------------------------------

.commodity.modelreport.priceBlack76:{[optionSetup;modelInputs]
    if[not `black76 in key modelInputs;
        :.commodity.modelreport.__failRow[`black76;"modelInputs missing key black76"]];
    inputs:modelInputs`black76;
    tryFn:{[setupDict;inputDict]
        if[not `volatility in key inputDict; '"modelInputs.black76 missing volatility"];
        volVal:inputDict`volatility;
        fwdVal:$[`forwardPrice in key inputDict; inputDict`forwardPrice; setupDict`forwardPrice];
        .commodity.black76.price[setupDict`optionType;fwdVal;setupDict`strikePrice;setupDict`expiry;volVal;setupDict`riskFreeRate]
        };
    priceVal:@[tryFn[optionSetup;];inputs;{x}];
    if[10h=type priceVal;
        :.commodity.modelreport.__failRow[`black76;priceVal]];
    .commodity.modelreport.__priceRow[`black76;priceVal;0Nf;0Nf;0Nf;`OK;""]
 };

.commodity.modelreport.priceSchwartz:{[optionSetup;modelInputs]
    if[not `schwartz in key modelInputs;
        :.commodity.modelreport.__failRow[`schwartz;"modelInputs missing key schwartz"]];
    inputs:modelInputs`schwartz;
    tryFn:{[setupDict;inputDict]
        if[not all `x0`params in key inputDict; '"modelInputs.schwartz missing x0 or params"];
        .commodity.schwartz.europeanOptionPrice[setupDict`optionType;inputDict`x0;setupDict`strikePrice;setupDict`expiry;setupDict`riskFreeRate;inputDict`params]
        };
    priceVal:@[tryFn[optionSetup;];inputs;{x}];
    if[10h=type priceVal;
        :.commodity.modelreport.__failRow[`schwartz;priceVal]];
    .commodity.modelreport.__priceRow[`schwartz;priceVal;0Nf;0Nf;0Nf;`OK;""]
 };

.commodity.modelreport.priceSchwartz2:{[optionSetup;modelInputs]
    if[not `schwartz2 in key modelInputs;
        :.commodity.modelreport.__failRow[`schwartz2;"modelInputs missing key schwartz2"]];
    inputs:modelInputs`schwartz2;
    tryFn:{[setupDict;inputDict]
        if[not all `shortFactor0`longFactor0`params in key inputDict;
            '"modelInputs.schwartz2 missing shortFactor0, longFactor0, or params"];
        .commodity.schwartz2.europeanOptionPrice[setupDict`optionType;inputDict`shortFactor0;inputDict`longFactor0;setupDict`strikePrice;setupDict`expiry;setupDict`riskFreeRate;inputDict`params]
        };
    priceVal:@[tryFn[optionSetup;];inputs;{x}];
    if[10h=type priceVal;
        :.commodity.modelreport.__failRow[`schwartz2;priceVal]];
    .commodity.modelreport.__priceRow[`schwartz2;priceVal;0Nf;0Nf;0Nf;`OK;""]
 };

.commodity.modelreport.priceMrJump:{[optionSetup;modelInputs]
    if[not `mrjump in key modelInputs;
        :.commodity.modelreport.__failRow[`mrjump;"modelInputs missing key mrjump"]];
    inputs:modelInputs`mrjump;
    tryFn:{[setupDict;inputDict]
        if[not all `x0`params`mcConfig in key inputDict;
            '"modelInputs.mrjump missing x0, params, or mcConfig"];
        .commodity.mrjump.europeanOptionPriceMC[setupDict`optionType;inputDict`x0;setupDict`strikePrice;setupDict`expiry;setupDict`riskFreeRate;inputDict`params;inputDict`mcConfig]
        };
    mcDict:@[tryFn[optionSetup;];inputs;{x}];
    if[10h=type mcDict;
        :.commodity.modelreport.__failRow[`mrjump;mcDict]];
    .commodity.modelreport.__priceRow[`mrjump;mcDict`price;mcDict`standardError;mcDict`lowerConfidence;mcDict`upperConfidence;`OK;""]
 };

/ --- Composite report API ------------------------------------------

.commodity.modelreport.priceAllModels:{[optionSetup;modelInputs]
    .commodity.modelreport.validateOptionSetup optionSetup;
    rowList:(
        .commodity.modelreport.priceBlack76[optionSetup;modelInputs];
        .commodity.modelreport.priceSchwartz[optionSetup;modelInputs];
        .commodity.modelreport.priceSchwartz2[optionSetup;modelInputs];
        .commodity.modelreport.priceMrJump[optionSetup;modelInputs]);
    .commodity.modelreport.__dictsToTable rowList
 };

.commodity.modelreport.comparisonSummary:{[priceTable]
    okMask:(priceTable`pricingStatus)=`OK;
    okRows:priceTable where okMask;
    modelCountValue:count priceTable;
    okModelCountValue:count okRows;
    baselineModel:`black76;
    baselinePriceVal:0Nf;
    baselineOkMask:okMask&(priceTable`modelName)=baselineModel;
    if[any baselineOkMask;
        baselineIdx:first where baselineOkMask;
        baselinePriceVal:priceTable[baselineIdx;`modelPrice]];
    avgPrice:0Nf; minPrice:0Nf; maxPrice:0Nf; priceRangeVal:0Nf;
    statusVal:`OK; errMsgVal:"";
    if[okModelCountValue>0;
        okPrices:okRows`modelPrice;
        avgPrice:avg okPrices;
        minPrice:min okPrices;
        maxPrice:max okPrices;
        priceRangeVal:maxPrice-minPrice];
    if[okModelCountValue=0;
        statusVal:`ERROR;
        errMsgVal:"All models failed pricing"];
    `baselineModel`baselinePrice`averageModelPrice`modelCount`okModelCount`minPrice`maxPrice`priceRange`status`errorMessage!(
        baselineModel;baselinePriceVal;avgPrice;modelCountValue;okModelCountValue;minPrice;maxPrice;priceRangeVal;statusVal;errMsgVal)
 };

.commodity.modelreport.priceDifferences:{[priceTable;baselineModel]
    summary:.commodity.modelreport.comparisonSummary priceTable;
    baselinePriceVal:summary`baselinePrice;
    avgPriceVal:summary`averageModelPrice;
    rowFn:{[row;baselineModelArg;baselinePriceArg;avgPriceArg]
        modelPriceVal:row`modelPrice;
        diffBaseline:modelPriceVal-baselinePriceArg;
        pctBaseline:$[(0f=baselinePriceArg)|null baselinePriceArg;
            0Nf;
            100f*diffBaseline%baselinePriceArg];
        diffAvg:modelPriceVal-avgPriceArg;
        `modelName`modelPrice`baselineModel`baselinePrice`differenceToBaseline`pctDifferenceToBaseline`differenceToAverage`pricingStatus`pricingErrorMessage!(
            row`modelName;modelPriceVal;baselineModelArg;baselinePriceArg;diffBaseline;pctBaseline;diffAvg;row`pricingStatus;row`pricingErrorMessage)
        };
    diffRows:rowFn[;baselineModel;baselinePriceVal;avgPriceVal] each priceTable;
    .commodity.modelreport.__dictsToTable diffRows
 };

/ --- Scenario application -----------------------------------------

.commodity.modelreport.scenarioShift:{[optionSetup;modelInputs;scenarioConfig]
    bumpedSetup:optionSetup;
    bumpedInputs:modelInputs;
    cfgKeys:key scenarioConfig;

    if[`forwardShiftPct in cfgKeys;
        scaleFn:.commodity.modelreport.__makeScaleFn scenarioConfig`forwardShiftPct;
        bumpedSetup:@[bumpedSetup;`forwardPrice;:;scaleFn bumpedSetup`forwardPrice]];

    if[`spotShiftPct in cfgKeys;
        shiftPctVal:scenarioConfig`spotShiftPct;
        scaleFn:.commodity.modelreport.__makeScaleFn shiftPctVal;
        addLogBumpFn:.commodity.modelreport.__makeAddFn log 1f+shiftPctVal;
        bumpedSetup:@[bumpedSetup;`spotPrice;:;scaleFn bumpedSetup`spotPrice];
        bumpedInputs:.commodity.modelreport.__bumpModelField[bumpedInputs;`schwartz;`x0;addLogBumpFn];
        bumpedInputs:.commodity.modelreport.__bumpModelField[bumpedInputs;`schwartz2;`longFactor0;addLogBumpFn];
        bumpedInputs:.commodity.modelreport.__bumpModelField[bumpedInputs;`mrjump;`x0;addLogBumpFn]];

    if[`volShift in cfgKeys;
        volBumpFn:.commodity.modelreport.__makeAddFn scenarioConfig`volShift;
        bumpedInputs:.commodity.modelreport.__bumpModelField[bumpedInputs;`black76;`volatility;volBumpFn];
        bumpedInputs:.commodity.modelreport.__bumpModelSubField[bumpedInputs;`schwartz;`params;`volatility;volBumpFn];
        bumpedInputs:.commodity.modelreport.__bumpModelSubField[bumpedInputs;`schwartz2;`params;`shortVolatility;volBumpFn];
        bumpedInputs:.commodity.modelreport.__bumpModelSubField[bumpedInputs;`mrjump;`params;`volatility;volBumpFn]];

    if[`longRunLogMeanShift in cfgKeys;
        meanShiftFn:.commodity.modelreport.__makeAddFn scenarioConfig`longRunLogMeanShift;
        bumpedInputs:.commodity.modelreport.__bumpModelSubField[bumpedInputs;`schwartz;`params;`longRunLogMean;meanShiftFn];
        bumpedInputs:.commodity.modelreport.__bumpModelField[bumpedInputs;`schwartz2;`longFactor0;meanShiftFn];
        bumpedInputs:.commodity.modelreport.__bumpModelSubField[bumpedInputs;`mrjump;`params;`longRunLogMean;meanShiftFn]];

    if[`jumpIntensityMultiplier in cfgKeys;
        mulFn:.commodity.modelreport.__makeMulFn scenarioConfig`jumpIntensityMultiplier;
        bumpedInputs:.commodity.modelreport.__bumpModelSubField[bumpedInputs;`mrjump;`params;`jumpIntensity;mulFn]];

    if[`jumpMeanShift in cfgKeys;
        jumpMeanFn:.commodity.modelreport.__makeAddFn scenarioConfig`jumpMeanShift;
        bumpedInputs:.commodity.modelreport.__bumpModelSubField[bumpedInputs;`mrjump;`params;`jumpMean;jumpMeanFn]];

    .commodity.modelreport.priceAllModels[bumpedSetup;bumpedInputs]
 };

/ --- Scenario PnL --------------------------------------------------

.commodity.modelreport.scenarioPnL:{[basePriceTable;scenarioPriceTable;quantity;contractMultiplier]
    rowFn:{[idx;baseTbl;scenTbl;qty;mult]
        baseRow:baseTbl idx;
        modelNameVal:baseRow`modelName;
        scenModels:scenTbl`modelName;
        scenMatch:scenModels?modelNameVal;
        scenarioPriceVal:0Nf;
        scenarioPnLVal:0Nf;
        rowStatus:`OK;
        errMsg:"";
        if[scenMatch>=count scenModels;
            rowStatus:`ERROR;
            errMsg:"No scenario row for model: ",string modelNameVal];
        if[scenMatch<count scenModels;
            scenRow:scenTbl scenMatch;
            scenarioPriceVal:scenRow`modelPrice;
            if[baseRow[`pricingStatus]=`ERROR;
                rowStatus:`ERROR;
                errMsg:"Base pricing failed"];
            if[scenRow[`pricingStatus]=`ERROR;
                rowStatus:`ERROR;
                errMsg:"Scenario pricing failed"];
            if[rowStatus=`OK;
                scenarioPnLVal:mult*qty*(scenarioPriceVal-baseRow`modelPrice)]];
        `modelName`basePrice`scenarioPrice`scenarioPnL`quantity`contractMultiplier`status`errorMessage!(
            modelNameVal;baseRow`modelPrice;scenarioPriceVal;scenarioPnLVal;qty;mult;rowStatus;errMsg)
        };
    rowDicts:rowFn[;basePriceTable;scenarioPriceTable;quantity;contractMultiplier] each til count basePriceTable;
    .commodity.modelreport.__dictsToTable rowDicts
 };

/ --- Convenience wrapper -------------------------------------------

.commodity.modelreport.runComparison:{[optionSetup;modelInputs;scenarioConfig;quantity;contractMultiplier]
    basePriceTableVal:.commodity.modelreport.priceAllModels[optionSetup;modelInputs];
    diffsVal:.commodity.modelreport.priceDifferences[basePriceTableVal;`black76];
    summaryVal:.commodity.modelreport.comparisonSummary basePriceTableVal;
    scenarioPriceTableVal:.commodity.modelreport.scenarioShift[optionSetup;modelInputs;scenarioConfig];
    pnLTableVal:.commodity.modelreport.scenarioPnL[basePriceTableVal;scenarioPriceTableVal;quantity;contractMultiplier];
    `basePrices`priceDifferences`summary`scenarioPrices`scenarioPnL!(
        basePriceTableVal;diffsVal;summaryVal;scenarioPriceTableVal;pnLTableVal)
 };

/ --- Optional ranking ---------------------------------------------

.commodity.modelreport.modelRanking:{[priceTable]
    okMask:(priceTable`pricingStatus)=`OK;
    okRows:priceTable where okMask;
    if[0=count okRows; :()];
    sortedRows:`modelPrice xdesc okRows;
    rankCount:count sortedRows;
    rankFn:{[idx;tbl]
        rowRecord:tbl idx;
        `rank`modelName`modelPrice!(1+idx;rowRecord`modelName;rowRecord`modelPrice)
        };
    rankDicts:rankFn[;sortedRows] each til rankCount;
    .commodity.modelreport.__dictsToTable rankDicts
 };

/ --- Greeks / sensitivities (v0.37) --------------------------------
/ Finite-difference sensitivities for each commodity model under a common
/ greekConfig. All four models use FD on their respective pricers, which keeps
/ the derivative units consistent (every column is dPrice / d(parameter)).
/ For mrjump the same mcConfig (and hence the same randomSeed) is reused
/ across base/up/down evaluations to apply common random numbers and reduce
/ FD noise.

.commodity.modelreport.defaultGreekConfig:{[]
    `spotBumpPct`volBump`meanReversionBump`longRunLogMeanBump`jumpIntensityBump`jumpMeanBump`correlationBump`useCentralDifference!(
        0.01;
        0.01;
        0.01;
        0.01;
        0.10;
        0.01;
        0.01;
        1b)
 };

.commodity.modelreport.bumpDict:{[paramDict;paramName;bumpAmount]
    if[not paramName in key paramDict; '"bumpDict missing key: ",string paramName];
    @[paramDict;paramName;+;bumpAmount]
 };

.commodity.modelreport.finiteDifference:{[basePrice;upPrice;downPrice;bumpAmount;useCentralDifference]
    if[0f=bumpAmount; '"finiteDifference bumpAmount must be non-zero"];
    $[useCentralDifference;
        (upPrice-downPrice)%(2f*bumpAmount);
        (upPrice-basePrice)%bumpAmount]
 };

/ --- Per-model greek wrappers ---

.commodity.modelreport.greeksBlack76:{[optionSetup;modelInputs;greekConfig]
    failRowDict:`modelName`basePrice`forwardDelta`volatilityVega`pricingStatus`pricingErrorMessage!(
        `black76;0Nf;0Nf;0Nf;`ERROR;"");
    bodyFn:{[setup;mi;cfg]
        if[not `black76 in key mi; '"modelInputs missing key black76"];
        inputs:mi`black76;
        if[not `volatility in key inputs; '"modelInputs.black76 missing volatility"];
        optType:setup`optionType;
        strikeVal:setup`strikePrice;
        expiryVal:setup`expiry;
        rateVal:setup`riskFreeRate;
        fwdVal:$[`forwardPrice in key inputs; inputs`forwardPrice; setup`forwardPrice];
        volVal:inputs`volatility;
        spotBumpPct:cfg`spotBumpPct;
        volBump:cfg`volBump;
        useCentral:cfg`useCentralDifference;
        fwdBumpAmount:fwdVal*spotBumpPct;
        basePrice:.commodity.black76.price[optType;fwdVal;strikeVal;expiryVal;volVal;rateVal];
        fwdUpPrice:.commodity.black76.price[optType;fwdVal+fwdBumpAmount;strikeVal;expiryVal;volVal;rateVal];
        fwdDownPrice:.commodity.black76.price[optType;fwdVal-fwdBumpAmount;strikeVal;expiryVal;volVal;rateVal];
        volUpPrice:.commodity.black76.price[optType;fwdVal;strikeVal;expiryVal;volVal+volBump;rateVal];
        volDownPrice:.commodity.black76.price[optType;fwdVal;strikeVal;expiryVal;volVal-volBump;rateVal];
        forwardDelta:.commodity.modelreport.finiteDifference[basePrice;fwdUpPrice;fwdDownPrice;fwdBumpAmount;useCentral];
        volatilityVega:.commodity.modelreport.finiteDifference[basePrice;volUpPrice;volDownPrice;volBump;useCentral];
        `basePrice`forwardDelta`volatilityVega!(basePrice;forwardDelta;volatilityVega)
        };
    result:.[bodyFn;(optionSetup;modelInputs;greekConfig);{x}];
    if[10h=type result;
        :.commodity.modelreport.__dictsToTable enlist @[failRowDict;`pricingErrorMessage;:;result]];
    okRow:`modelName`basePrice`forwardDelta`volatilityVega`pricingStatus`pricingErrorMessage!(
        `black76;result`basePrice;result`forwardDelta;result`volatilityVega;`OK;"");
    .commodity.modelreport.__dictsToTable enlist okRow
 };

.commodity.modelreport.greeksSchwartz:{[optionSetup;modelInputs;greekConfig]
    failRowDict:`modelName`basePrice`logStateDelta`volatilityVega`meanReversionSensitivity`longRunLogMeanSensitivity`pricingStatus`pricingErrorMessage!(
        `schwartz;0Nf;0Nf;0Nf;0Nf;0Nf;`ERROR;"");
    bodyFn:{[setup;mi;cfg]
        if[not `schwartz in key mi; '"modelInputs missing key schwartz"];
        inputs:mi`schwartz;
        if[not all `x0`params in key inputs; '"modelInputs.schwartz missing x0 or params"];
        optType:setup`optionType;
        strikeVal:setup`strikePrice;
        expiryVal:setup`expiry;
        rateVal:setup`riskFreeRate;
        x0:inputs`x0;
        params:inputs`params;
        spotBumpPct:cfg`spotBumpPct;
        logBump:log 1f+spotBumpPct;
        volBump:cfg`volBump;
        kappaBump:cfg`meanReversionBump;
        thetaBump:cfg`longRunLogMeanBump;
        useCentral:cfg`useCentralDifference;
        priceFn:.commodity.schwartz.europeanOptionPrice[optType;;strikeVal;expiryVal;rateVal;];
        basePrice:priceFn[x0;params];
        xUpPrice:priceFn[x0+logBump;params];
        xDownPrice:priceFn[x0-logBump;params];
        logStateDelta:.commodity.modelreport.finiteDifference[basePrice;xUpPrice;xDownPrice;logBump;useCentral];
        volUpPrice:priceFn[x0;.commodity.modelreport.bumpDict[params;`volatility;volBump]];
        volDownPrice:priceFn[x0;.commodity.modelreport.bumpDict[params;`volatility;neg volBump]];
        volatilityVega:.commodity.modelreport.finiteDifference[basePrice;volUpPrice;volDownPrice;volBump;useCentral];
        kappaUpPrice:priceFn[x0;.commodity.modelreport.bumpDict[params;`meanReversionSpeed;kappaBump]];
        kappaDownPrice:priceFn[x0;.commodity.modelreport.bumpDict[params;`meanReversionSpeed;neg kappaBump]];
        meanReversionSensitivity:.commodity.modelreport.finiteDifference[basePrice;kappaUpPrice;kappaDownPrice;kappaBump;useCentral];
        thetaUpPrice:priceFn[x0;.commodity.modelreport.bumpDict[params;`longRunLogMean;thetaBump]];
        thetaDownPrice:priceFn[x0;.commodity.modelreport.bumpDict[params;`longRunLogMean;neg thetaBump]];
        longRunLogMeanSensitivity:.commodity.modelreport.finiteDifference[basePrice;thetaUpPrice;thetaDownPrice;thetaBump;useCentral];
        `basePrice`logStateDelta`volatilityVega`meanReversionSensitivity`longRunLogMeanSensitivity!(
            basePrice;logStateDelta;volatilityVega;meanReversionSensitivity;longRunLogMeanSensitivity)
        };
    result:.[bodyFn;(optionSetup;modelInputs;greekConfig);{x}];
    if[10h=type result;
        :.commodity.modelreport.__dictsToTable enlist @[failRowDict;`pricingErrorMessage;:;result]];
    okRow:`modelName`basePrice`logStateDelta`volatilityVega`meanReversionSensitivity`longRunLogMeanSensitivity`pricingStatus`pricingErrorMessage!(
        `schwartz;result`basePrice;result`logStateDelta;result`volatilityVega;result`meanReversionSensitivity;result`longRunLogMeanSensitivity;`OK;"");
    .commodity.modelreport.__dictsToTable enlist okRow
 };

.commodity.modelreport.greeksSchwartz2:{[optionSetup;modelInputs;greekConfig]
    failRowDict:`modelName`basePrice`shortFactorSensitivity`longFactorSensitivity`shortVolSensitivity`longVolSensitivity`meanReversionSensitivity`correlationSensitivity`pricingStatus`pricingErrorMessage!(
        `schwartz2;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;`ERROR;"");
    bodyFn:{[setup;mi;cfg]
        if[not `schwartz2 in key mi; '"modelInputs missing key schwartz2"];
        inputs:mi`schwartz2;
        if[not all `shortFactor0`longFactor0`params in key inputs;
            '"modelInputs.schwartz2 missing shortFactor0, longFactor0, or params"];
        optType:setup`optionType;
        strikeVal:setup`strikePrice;
        expiryVal:setup`expiry;
        rateVal:setup`riskFreeRate;
        shortFactor0:inputs`shortFactor0;
        longFactor0:inputs`longFactor0;
        params:inputs`params;
        spotBumpPct:cfg`spotBumpPct;
        logBump:log 1f+spotBumpPct;
        volBump:cfg`volBump;
        kappaBump:cfg`meanReversionBump;
        rhoBump:cfg`correlationBump;
        useCentral:cfg`useCentralDifference;
        rhoVal:params`correlation;
        if[(rhoVal+rhoBump)>=1f; '"schwartz2 correlation bump would exceed +1; reduce correlationBump"];
        if[(rhoVal-rhoBump)<=-1f; '"schwartz2 correlation bump would breach -1; reduce correlationBump"];
        priceFn:.commodity.schwartz2.europeanOptionPrice[optType;;;strikeVal;expiryVal;rateVal;];
        basePrice:priceFn[shortFactor0;longFactor0;params];
        shortUpPrice:priceFn[shortFactor0+logBump;longFactor0;params];
        shortDownPrice:priceFn[shortFactor0-logBump;longFactor0;params];
        shortFactorSensitivity:.commodity.modelreport.finiteDifference[basePrice;shortUpPrice;shortDownPrice;logBump;useCentral];
        longUpPrice:priceFn[shortFactor0;longFactor0+logBump;params];
        longDownPrice:priceFn[shortFactor0;longFactor0-logBump;params];
        longFactorSensitivity:.commodity.modelreport.finiteDifference[basePrice;longUpPrice;longDownPrice;logBump;useCentral];
        shortVolUpPrice:priceFn[shortFactor0;longFactor0;.commodity.modelreport.bumpDict[params;`shortVolatility;volBump]];
        shortVolDownPrice:priceFn[shortFactor0;longFactor0;.commodity.modelreport.bumpDict[params;`shortVolatility;neg volBump]];
        shortVolSensitivity:.commodity.modelreport.finiteDifference[basePrice;shortVolUpPrice;shortVolDownPrice;volBump;useCentral];
        longVolUpPrice:priceFn[shortFactor0;longFactor0;.commodity.modelreport.bumpDict[params;`longVolatility;volBump]];
        longVolDownPrice:priceFn[shortFactor0;longFactor0;.commodity.modelreport.bumpDict[params;`longVolatility;neg volBump]];
        longVolSensitivity:.commodity.modelreport.finiteDifference[basePrice;longVolUpPrice;longVolDownPrice;volBump;useCentral];
        kappaUpPrice:priceFn[shortFactor0;longFactor0;.commodity.modelreport.bumpDict[params;`meanReversionSpeed;kappaBump]];
        kappaDownPrice:priceFn[shortFactor0;longFactor0;.commodity.modelreport.bumpDict[params;`meanReversionSpeed;neg kappaBump]];
        meanReversionSensitivity:.commodity.modelreport.finiteDifference[basePrice;kappaUpPrice;kappaDownPrice;kappaBump;useCentral];
        rhoUpPrice:priceFn[shortFactor0;longFactor0;.commodity.modelreport.bumpDict[params;`correlation;rhoBump]];
        rhoDownPrice:priceFn[shortFactor0;longFactor0;.commodity.modelreport.bumpDict[params;`correlation;neg rhoBump]];
        correlationSensitivity:.commodity.modelreport.finiteDifference[basePrice;rhoUpPrice;rhoDownPrice;rhoBump;useCentral];
        `basePrice`shortFactorSensitivity`longFactorSensitivity`shortVolSensitivity`longVolSensitivity`meanReversionSensitivity`correlationSensitivity!(
            basePrice;shortFactorSensitivity;longFactorSensitivity;shortVolSensitivity;longVolSensitivity;meanReversionSensitivity;correlationSensitivity)
        };
    result:.[bodyFn;(optionSetup;modelInputs;greekConfig);{x}];
    if[10h=type result;
        :.commodity.modelreport.__dictsToTable enlist @[failRowDict;`pricingErrorMessage;:;result]];
    okRow:`modelName`basePrice`shortFactorSensitivity`longFactorSensitivity`shortVolSensitivity`longVolSensitivity`meanReversionSensitivity`correlationSensitivity`pricingStatus`pricingErrorMessage!(
        `schwartz2;result`basePrice;result`shortFactorSensitivity;result`longFactorSensitivity;result`shortVolSensitivity;result`longVolSensitivity;result`meanReversionSensitivity;result`correlationSensitivity;`OK;"");
    .commodity.modelreport.__dictsToTable enlist okRow
 };

.commodity.modelreport.greeksMrJump:{[optionSetup;modelInputs;greekConfig]
    failRowDict:`modelName`basePrice`logStateDelta`volatilityVega`meanReversionSensitivity`longRunLogMeanSensitivity`jumpIntensitySensitivity`jumpMeanSensitivity`standardError`pricingStatus`pricingErrorMessage!(
        `mrjump;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;`ERROR;"");
    bodyFn:{[setup;mi;cfg]
        if[not `mrjump in key mi; '"modelInputs missing key mrjump"];
        inputs:mi`mrjump;
        if[not all `x0`params`mcConfig in key inputs;
            '"modelInputs.mrjump missing x0, params, or mcConfig"];
        optType:setup`optionType;
        strikeVal:setup`strikePrice;
        expiryVal:setup`expiry;
        rateVal:setup`riskFreeRate;
        x0:inputs`x0;
        params:inputs`params;
        mcCfg:inputs`mcConfig;
        spotBumpPct:cfg`spotBumpPct;
        logBump:log 1f+spotBumpPct;
        volBump:cfg`volBump;
        kappaBump:cfg`meanReversionBump;
        thetaBump:cfg`longRunLogMeanBump;
        lambdaBump:cfg`jumpIntensityBump;
        jumpMeanBumpVal:cfg`jumpMeanBump;
        useCentral:cfg`useCentralDifference;
        if[(params`jumpIntensity)<lambdaBump;
            '"mrjump greeks: jumpIntensityBump would push jumpIntensity below zero; choose a smaller bump"];
        priceFn:{[ot;k;t;r;mc;xVal;pVal]
            (.commodity.mrjump.europeanOptionPriceMC[ot;xVal;k;t;r;pVal;mc])`price
            }[optType;strikeVal;expiryVal;rateVal;mcCfg;;];
        baseMcResult:.commodity.mrjump.europeanOptionPriceMC[optType;x0;strikeVal;expiryVal;rateVal;params;mcCfg];
        basePrice:baseMcResult`price;
        baseSE:baseMcResult`standardError;
        xUpPrice:priceFn[x0+logBump;params];
        xDownPrice:priceFn[x0-logBump;params];
        logStateDelta:.commodity.modelreport.finiteDifference[basePrice;xUpPrice;xDownPrice;logBump;useCentral];
        volUpPrice:priceFn[x0;.commodity.modelreport.bumpDict[params;`volatility;volBump]];
        volDownPrice:priceFn[x0;.commodity.modelreport.bumpDict[params;`volatility;neg volBump]];
        volatilityVega:.commodity.modelreport.finiteDifference[basePrice;volUpPrice;volDownPrice;volBump;useCentral];
        kappaUpPrice:priceFn[x0;.commodity.modelreport.bumpDict[params;`meanReversionSpeed;kappaBump]];
        kappaDownPrice:priceFn[x0;.commodity.modelreport.bumpDict[params;`meanReversionSpeed;neg kappaBump]];
        meanReversionSensitivity:.commodity.modelreport.finiteDifference[basePrice;kappaUpPrice;kappaDownPrice;kappaBump;useCentral];
        thetaUpPrice:priceFn[x0;.commodity.modelreport.bumpDict[params;`longRunLogMean;thetaBump]];
        thetaDownPrice:priceFn[x0;.commodity.modelreport.bumpDict[params;`longRunLogMean;neg thetaBump]];
        longRunLogMeanSensitivity:.commodity.modelreport.finiteDifference[basePrice;thetaUpPrice;thetaDownPrice;thetaBump;useCentral];
        lambdaUpPrice:priceFn[x0;.commodity.modelreport.bumpDict[params;`jumpIntensity;lambdaBump]];
        lambdaDownPrice:priceFn[x0;.commodity.modelreport.bumpDict[params;`jumpIntensity;neg lambdaBump]];
        jumpIntensitySensitivity:.commodity.modelreport.finiteDifference[basePrice;lambdaUpPrice;lambdaDownPrice;lambdaBump;useCentral];
        jmUpPrice:priceFn[x0;.commodity.modelreport.bumpDict[params;`jumpMean;jumpMeanBumpVal]];
        jmDownPrice:priceFn[x0;.commodity.modelreport.bumpDict[params;`jumpMean;neg jumpMeanBumpVal]];
        jumpMeanSensitivity:.commodity.modelreport.finiteDifference[basePrice;jmUpPrice;jmDownPrice;jumpMeanBumpVal;useCentral];
        `basePrice`logStateDelta`volatilityVega`meanReversionSensitivity`longRunLogMeanSensitivity`jumpIntensitySensitivity`jumpMeanSensitivity`standardError!(
            basePrice;logStateDelta;volatilityVega;meanReversionSensitivity;longRunLogMeanSensitivity;jumpIntensitySensitivity;jumpMeanSensitivity;baseSE)
        };
    result:.[bodyFn;(optionSetup;modelInputs;greekConfig);{x}];
    if[10h=type result;
        :.commodity.modelreport.__dictsToTable enlist @[failRowDict;`pricingErrorMessage;:;result]];
    okRow:`modelName`basePrice`logStateDelta`volatilityVega`meanReversionSensitivity`longRunLogMeanSensitivity`jumpIntensitySensitivity`jumpMeanSensitivity`standardError`pricingStatus`pricingErrorMessage!(
        `mrjump;result`basePrice;result`logStateDelta;result`volatilityVega;result`meanReversionSensitivity;result`longRunLogMeanSensitivity;result`jumpIntensitySensitivity;result`jumpMeanSensitivity;result`standardError;`OK;"");
    .commodity.modelreport.__dictsToTable enlist okRow
 };

/ --- Composite greeks API ---

.commodity.modelreport.greeksAllModels:{[optionSetup;modelInputs;greekConfig]
    b76Tbl:.commodity.modelreport.greeksBlack76[optionSetup;modelInputs;greekConfig];
    schTbl:.commodity.modelreport.greeksSchwartz[optionSetup;modelInputs;greekConfig];
    sch2Tbl:.commodity.modelreport.greeksSchwartz2[optionSetup;modelInputs;greekConfig];
    mrjTbl:.commodity.modelreport.greeksMrJump[optionSetup;modelInputs;greekConfig];
    b76Tbl uj schTbl uj sch2Tbl uj mrjTbl
 };

.commodity.modelreport.greeksSummary:{[greeksTable]
    okMask:(greeksTable`pricingStatus)=`OK;
    okRows:greeksTable where okMask;
    modelCountValue:count greeksTable;
    okModelCountValue:count okRows;
    statusVal:`OK; errMsgVal:"";
    if[okModelCountValue=0;
        statusVal:`ERROR;
        errMsgVal:"All model greeks failed"];
    maxPrimaryDeltaValue:0Nf;
    maxVolatilityVegaValue:0Nf;
    maxJumpIntensitySensitivityValue:0Nf;
    if[okModelCountValue>0;
        availableCols:cols okRows;
        primaryDeltaFn:{[rowDict]
            modelNameVal:rowDict`modelName;
            $[modelNameVal=`black76; rowDict`forwardDelta;
              modelNameVal=`schwartz; rowDict`logStateDelta;
              modelNameVal=`schwartz2; rowDict`longFactorSensitivity;
              modelNameVal=`mrjump; rowDict`logStateDelta;
              0Nf]
            };
        primaryDeltas:primaryDeltaFn each okRows;
        nonNullDeltas:primaryDeltas where not null primaryDeltas;
        if[0<count nonNullDeltas; maxPrimaryDeltaValue:max abs nonNullDeltas];
        if[`volatilityVega in availableCols;
            volVegas:okRows`volatilityVega;
            nonNullVegas:volVegas where not null volVegas;
            if[0<count nonNullVegas; maxVolatilityVegaValue:max abs nonNullVegas]];
        if[`jumpIntensitySensitivity in availableCols;
            jumpSens:okRows`jumpIntensitySensitivity;
            nonNullJumpSens:jumpSens where not null jumpSens;
            if[0<count nonNullJumpSens; maxJumpIntensitySensitivityValue:max abs nonNullJumpSens]]];
    `modelCount`okModelCount`maxPrimaryDelta`maxVolatilityVega`maxJumpIntensitySensitivity`status`errorMessage!(
        modelCountValue;okModelCountValue;maxPrimaryDeltaValue;maxVolatilityVegaValue;maxJumpIntensitySensitivityValue;statusVal;errMsgVal)
 };

.commodity.modelreport.runComparisonWithGreeks:{[optionSetup;modelInputs;scenarioConfig;quantity;contractMultiplier;greekConfig]
    baseReport:.commodity.modelreport.runComparison[optionSetup;modelInputs;scenarioConfig;quantity;contractMultiplier];
    greeksTbl:.commodity.modelreport.greeksAllModels[optionSetup;modelInputs;greekConfig];
    greeksSum:.commodity.modelreport.greeksSummary greeksTbl;
    baseReport,`greeksTable`greeksSummary!(greeksTbl;greeksSum)
 };

/ --- Disagreement risk (v0.38) -------------------------------------
/ Quantifies how much price/risk depends on commodity model choice.
/ Public:  defaultDisagreementConfig, primarySensitivity, priceDisagreement,
/          greeksDisagreement, scenarioDisagreement, disagreementAlerts,
/          modelDisagreementReport, runComparisonRisk.
/ Private: __primaryFieldFor, __resolveStatus, __alertRow.
/ Each disagreement dict carries its own threshold values so the alert table
/ can be built without re-passing the disagreementConfig.

.commodity.modelreport.defaultDisagreementConfig:{[]
    `priceRangeAbsThreshold`priceRangePctThreshold`primaryDeltaRangeThreshold`volatilityVegaRangeThreshold`scenarioPnlRangeThreshold`jumpSensitivityThreshold`minimumOkModels!(
        5f;
        0.10;
        10f;
        10f;
        1000f;
        5f;
        2)
 };

.commodity.modelreport.__primaryFieldFor:{[modelName]
    $[modelName=`black76; `forwardDelta;
      modelName=`schwartz; `logStateDelta;
      modelName=`schwartz2; `longFactorSensitivity;
      modelName=`mrjump; `logStateDelta;
      `]
 };

/ Common status resolution: zero OK rows -> ERROR; some but below minimum -> warning; otherwise OK.
.commodity.modelreport.__resolveStatus:{[okCount;minimumRequired;errorContextLabel]
    if[okCount=0;
        :`status`errorMessage!(`ERROR;"No OK rows for ",errorContextLabel)];
    if[okCount<minimumRequired;
        :`status`errorMessage!(`warning;errorContextLabel," has ",string[okCount]," OK rows < minimumOkModels ",string minimumRequired)];
    `status`errorMessage!(`OK;"")
 };

.commodity.modelreport.__alertRow:{[alertType;alertFlag;metricValue;threshold;message]
    severity:$[alertFlag;`warning;`OK];
    `alertType`alertFlag`metricValue`threshold`severity`message!(
        alertType;alertFlag;metricValue;threshold;severity;message)
 };

/ --- Primary sensitivity extraction ---

.commodity.modelreport.primarySensitivity:{[greeksTable]
    rowFn:{[rowDict]
        modelNameVal:rowDict`modelName;
        primaryField:.commodity.modelreport.__primaryFieldFor modelNameVal;
        rowStatus:rowDict`pricingStatus;
        rowErrMsg:rowDict`pricingErrorMessage;
        primaryVal:0Nf;
        if[(rowStatus=`OK)&(not primaryField=`)&(primaryField in key rowDict);
            candidate:rowDict primaryField;
            if[not null candidate; primaryVal:candidate]];
        `modelName`primarySensitivityName`primarySensitivity`pricingStatus`pricingErrorMessage!(
            modelNameVal;primaryField;primaryVal;rowStatus;rowErrMsg)
        };
    rowDicts:rowFn each greeksTable;
    .commodity.modelreport.__dictsToTable rowDicts
 };

/ --- Price disagreement ---

.commodity.modelreport.priceDisagreement:{[priceTable;disagreementConfig]
    minimumRequired:disagreementConfig`minimumOkModels;
    absThr:disagreementConfig`priceRangeAbsThreshold;
    pctThr:disagreementConfig`priceRangePctThreshold;
    okMask:(priceTable`pricingStatus)=`OK;
    okRows:priceTable where okMask;
    okPrices:okRows`modelPrice;
    validMask:not null okPrices;
    validPrices:okPrices where validMask;
    okModelCountValue:count validPrices;

    minPriceVal:0Nf; maxPriceVal:0Nf; averagePriceVal:0Nf;
    priceRangeVal:0Nf; priceRangePctVal:0Nf;
    priceRangeAlertVal:0b;
    if[okModelCountValue>0;
        minPriceVal:min validPrices;
        maxPriceVal:max validPrices;
        averagePriceVal:avg validPrices;
        priceRangeVal:maxPriceVal-minPriceVal;
        priceRangePctVal:$[0f=averagePriceVal;0Nf;priceRangeVal%abs averagePriceVal];
        absAlert:priceRangeVal>absThr;
        pctAlert:(not null priceRangePctVal)&priceRangePctVal>pctThr;
        priceRangeAlertVal:absAlert|pctAlert];

    statusDict:.commodity.modelreport.__resolveStatus[okModelCountValue;minimumRequired;"priceDisagreement"];
    `okModelCount`minPrice`maxPrice`averagePrice`priceRange`priceRangePct`priceRangeAlert`priceRangeAbsThreshold`priceRangePctThreshold`status`errorMessage!(
        okModelCountValue;minPriceVal;maxPriceVal;averagePriceVal;priceRangeVal;priceRangePctVal;priceRangeAlertVal;absThr;pctThr;statusDict`status;statusDict`errorMessage)
 };

/ --- Greeks disagreement ---

.commodity.modelreport.greeksDisagreement:{[greeksTable;disagreementConfig]
    minimumRequired:disagreementConfig`minimumOkModels;
    primaryThr:disagreementConfig`primaryDeltaRangeThreshold;
    volVegaThr:disagreementConfig`volatilityVegaRangeThreshold;
    jumpThr:disagreementConfig`jumpSensitivityThreshold;

    primaryTbl:.commodity.modelreport.primarySensitivity greeksTable;
    primaryOkMask:(primaryTbl`pricingStatus)=`OK;
    primaryOkRows:primaryTbl where primaryOkMask;
    primaryVals:primaryOkRows`primarySensitivity;
    validPrimaryMask:not null primaryVals;
    validPrimary:primaryVals where validPrimaryMask;
    okModelCountValue:count validPrimary;

    primaryMinVal:0Nf; primaryMaxVal:0Nf; primaryRangeVal:0Nf;
    primaryAlertVal:0b;
    if[okModelCountValue>0;
        primaryMinVal:min validPrimary;
        primaryMaxVal:max validPrimary;
        primaryRangeVal:primaryMaxVal-primaryMinVal;
        primaryAlertVal:primaryRangeVal>primaryThr];

    greekOkMask:(greeksTable`pricingStatus)=`OK;
    greekOkRows:greeksTable where greekOkMask;
    availableGreekCols:cols greekOkRows;

    volVegaMinVal:0Nf; volVegaMaxVal:0Nf; volVegaRangeVal:0Nf;
    volVegaAlertVal:0b;
    if[`volatilityVega in availableGreekCols;
        volVegasFull:greekOkRows`volatilityVega;
        validVegaMask:not null volVegasFull;
        validVegas:volVegasFull where validVegaMask;
        if[0<count validVegas;
            volVegaMinVal:min validVegas;
            volVegaMaxVal:max validVegas;
            volVegaRangeVal:volVegaMaxVal-volVegaMinVal;
            volVegaAlertVal:volVegaRangeVal>volVegaThr]];

    jumpIntSensVal:0Nf;
    jumpSensAlertVal:0b;
    if[`jumpIntensitySensitivity in availableGreekCols;
        mrjumpMask:(greekOkRows`modelName)=`mrjump;
        if[any mrjumpMask;
            mrjumpRowSelect:greekOkRows first where mrjumpMask;
            candidate:mrjumpRowSelect`jumpIntensitySensitivity;
            if[not null candidate;
                jumpIntSensVal:candidate;
                jumpSensAlertVal:(abs jumpIntSensVal)>jumpThr]]];

    statusDict:.commodity.modelreport.__resolveStatus[okModelCountValue;minimumRequired;"greeksDisagreement"];
    `okModelCount`primarySensitivityMin`primarySensitivityMax`primarySensitivityRange`primarySensitivityAlert`volatilityVegaMin`volatilityVegaMax`volatilityVegaRange`volatilityVegaAlert`jumpIntensitySensitivity`jumpSensitivityAlert`primaryDeltaRangeThreshold`volatilityVegaRangeThreshold`jumpSensitivityThreshold`status`errorMessage!(
        okModelCountValue;primaryMinVal;primaryMaxVal;primaryRangeVal;primaryAlertVal;volVegaMinVal;volVegaMaxVal;volVegaRangeVal;volVegaAlertVal;jumpIntSensVal;jumpSensAlertVal;primaryThr;volVegaThr;jumpThr;statusDict`status;statusDict`errorMessage)
 };

/ --- Scenario PnL disagreement ---

.commodity.modelreport.scenarioDisagreement:{[scenarioPnlTable;disagreementConfig]
    minimumRequired:disagreementConfig`minimumOkModels;
    pnlThr:disagreementConfig`scenarioPnlRangeThreshold;
    okMask:(scenarioPnlTable`status)=`OK;
    okRows:scenarioPnlTable where okMask;
    pnLs:okRows`scenarioPnL;
    validMask:not null pnLs;
    validPnLs:pnLs where validMask;
    okModelCountValue:count validPnLs;

    minPnlVal:0Nf; maxPnlVal:0Nf; avgPnlVal:0Nf; rangeVal:0Nf;
    alertVal:0b;
    if[okModelCountValue>0;
        minPnlVal:min validPnLs;
        maxPnlVal:max validPnLs;
        avgPnlVal:avg validPnLs;
        rangeVal:maxPnlVal-minPnlVal;
        alertVal:rangeVal>pnlThr];

    statusDict:.commodity.modelreport.__resolveStatus[okModelCountValue;minimumRequired;"scenarioDisagreement"];
    `okModelCount`minScenarioPnl`maxScenarioPnl`averageScenarioPnl`scenarioPnlRange`scenarioPnlAlert`scenarioPnlRangeThreshold`status`errorMessage!(
        okModelCountValue;minPnlVal;maxPnlVal;avgPnlVal;rangeVal;alertVal;pnlThr;statusDict`status;statusDict`errorMessage)
 };

/ --- Alert composition ---

.commodity.modelreport.disagreementAlerts:{[priceDisagreement;greeksDisagreement;scenarioDisagreement]
    rowList:(
        .commodity.modelreport.__alertRow[`priceRangeAbs;
            priceDisagreement[`priceRange]>priceDisagreement`priceRangeAbsThreshold;
            priceDisagreement`priceRange;
            priceDisagreement`priceRangeAbsThreshold;
            "Absolute price range across models"];
        .commodity.modelreport.__alertRow[`priceRangePct;
            (not null priceDisagreement`priceRangePct)&priceDisagreement[`priceRangePct]>priceDisagreement`priceRangePctThreshold;
            priceDisagreement`priceRangePct;
            priceDisagreement`priceRangePctThreshold;
            "Relative price range across models"];
        .commodity.modelreport.__alertRow[`primarySensitivityRange;
            greeksDisagreement`primarySensitivityAlert;
            greeksDisagreement`primarySensitivityRange;
            greeksDisagreement`primaryDeltaRangeThreshold;
            "Primary-sensitivity range across models"];
        .commodity.modelreport.__alertRow[`volatilityVegaRange;
            greeksDisagreement`volatilityVegaAlert;
            greeksDisagreement`volatilityVegaRange;
            greeksDisagreement`volatilityVegaRangeThreshold;
            "Volatility-vega range across models"];
        .commodity.modelreport.__alertRow[`scenarioPnlRange;
            scenarioDisagreement`scenarioPnlAlert;
            scenarioDisagreement`scenarioPnlRange;
            scenarioDisagreement`scenarioPnlRangeThreshold;
            "Scenario PnL range across models"];
        .commodity.modelreport.__alertRow[`jumpSensitivity;
            greeksDisagreement`jumpSensitivityAlert;
            greeksDisagreement`jumpIntensitySensitivity;
            greeksDisagreement`jumpSensitivityThreshold;
            "mrjump jump-intensity sensitivity magnitude"]);
    .commodity.modelreport.__dictsToTable rowList
 };

/ --- Full report wrappers ---

.commodity.modelreport.modelDisagreementReport:{[comparisonResult;disagreementConfig]
    basePricesTable:comparisonResult`basePrices;
    greeksTbl:comparisonResult`greeksTable;
    scenarioPnLTbl:comparisonResult`scenarioPnL;
    priceDis:.commodity.modelreport.priceDisagreement[basePricesTable;disagreementConfig];
    greeksDis:.commodity.modelreport.greeksDisagreement[greeksTbl;disagreementConfig];
    scenarioDis:.commodity.modelreport.scenarioDisagreement[scenarioPnLTbl;disagreementConfig];
    alertsTbl:.commodity.modelreport.disagreementAlerts[priceDis;greeksDis;scenarioDis];
    `priceDisagreement`greeksDisagreement`scenarioDisagreement`alerts!(
        priceDis;greeksDis;scenarioDis;alertsTbl)
 };

.commodity.modelreport.runComparisonRisk:{[optionSetup;modelInputs;scenarioConfig;quantity;contractMultiplier;greekConfig;disagreementConfig]
    comparisonResult:.commodity.modelreport.runComparisonWithGreeks[optionSetup;modelInputs;scenarioConfig;quantity;contractMultiplier;greekConfig];
    disagreementResult:.commodity.modelreport.modelDisagreementReport[comparisonResult;disagreementConfig];
    `comparison`disagreement!(comparisonResult;disagreementResult)
 };
