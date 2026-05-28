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
