/ modelcheck.q - model limit and consistency checks (v0.26)

/ --- Default configs ---

.modelcheck.__defaultTrade:{[]
    `tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
        1;`TEST;`equityOption;`european;`call;100f;1f;1f;`none;0Nf;0f)
 };

.modelcheck.__defaultMarket:{[]
    `underlying`spot`riskFreeRate`dividendYield`volatility!(`TEST;100f;0.05;0f;0.2)
 };

.modelcheck.__defaultMcConfig:{[]
    `pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(10000;50;42;0b;0b;0.95)
 };


.modelcheck.__getOr:{[dict;keyName;defaultVal]
    if[not 99h=type dict; :defaultVal];
    if[0=count dict; :defaultVal];
    if[keyName in key dict; :dict keyName];
    defaultVal
 };

/ --- Heston -> BS ---

.modelcheck.checkHestonBlackScholesLimit:{[trade;marketData;configDict]
    checkName:"HestonToBS";
    mcConfig:.modelcheck.__getOr[configDict;`mcConfig;.modelcheck.__defaultMcConfig[]];
    absTol:.modelcheck.__getOr[configDict;`absoluteTolerance;0.5];
    relTol:.modelcheck.__getOr[configDict;`relativeTolerance;0.05];
    bsPrice:.validation.blackScholesClosedForm[trade`optionType;marketData`spot;trade`strike;trade`expiry;marketData`riskFreeRate;marketData`dividendYield;marketData`volatility];
    sigmaSquared:marketData[`volatility]*marketData`volatility;
    hestonParams:`initialVariance`longRunVariance`meanReversion`volOfVol`correlation`riskFreeRate`dividendYield!(
        sigmaSquared;sigmaSquared;2.0;0.0;0.0;marketData`riskFreeRate;marketData`dividendYield);
    hestonConfig:`mcConfig`hestonParams`modelType!(mcConfig;hestonParams;`heston);
    hestonFn:{.heston.priceEuropean[x 0;x 1;x 2]};
    hestonResult:@[hestonFn;(trade;marketData;hestonConfig);{x}];
    if[10h=type hestonResult; :.limitcheck.errorRow[checkName;`heston;`blackScholes;hestonResult]];
    .limitcheck.resultRow[checkName;`heston;`blackScholes;hestonResult`unitPrice;bsPrice;absTol;relTol]
 };

/ --- Merton -> BS ---

.modelcheck.checkMertonBlackScholesLimit:{[trade;marketData;configDict]
    checkName:"MertonToBS";
    absTol:.modelcheck.__getOr[configDict;`absoluteTolerance;0.01];
    relTol:.modelcheck.__getOr[configDict;`relativeTolerance;0.001];
    bsPrice:.validation.blackScholesClosedForm[trade`optionType;marketData`spot;trade`strike;trade`expiry;marketData`riskFreeRate;marketData`dividendYield;marketData`volatility];
    mertonParams:`volatility`jumpIntensity`jumpMean`jumpVolatility`riskFreeRate`dividendYield!(
        marketData`volatility;0.0;0.0;0.0;marketData`riskFreeRate;marketData`dividendYield);
    mertonPrice:.merton.priceEuropeanSeries[trade`optionType;marketData`spot;trade`strike;trade`expiry;mertonParams;30];
    .limitcheck.resultRow[checkName;`merton;`blackScholes;mertonPrice;bsPrice;absTol;relTol]
 };

/ --- Bates -> Heston ---

.modelcheck.checkBatesHestonLimit:{[trade;marketData;configDict]
    checkName:"BatesToHeston";
    mcConfig:.modelcheck.__getOr[configDict;`mcConfig;.modelcheck.__defaultMcConfig[]];
    absTol:.modelcheck.__getOr[configDict;`absoluteTolerance;1.0];
    relTol:.modelcheck.__getOr[configDict;`relativeTolerance;0.1];
    hestonParams:`initialVariance`longRunVariance`meanReversion`volOfVol`correlation`riskFreeRate`dividendYield!(
        0.04;0.04;2.0;0.3;-0.7;marketData`riskFreeRate;marketData`dividendYield);
    hestonConfig:`mcConfig`hestonParams`modelType!(mcConfig;hestonParams;`heston);
    hestonFn:{.heston.priceEuropean[x 0;x 1;x 2]};
    hestonResult:@[hestonFn;(trade;marketData;hestonConfig);{x}];
    if[10h=type hestonResult; :.limitcheck.errorRow[checkName;`bates;`heston;hestonResult]];
    batesParams:hestonParams,`jumpIntensity`jumpMean`jumpVolatility!(0.0;0.0;0.0);
    batesConfig:`mcConfig`batesParams!(mcConfig;batesParams);
    batesFn:{.bates.priceEuropean[x 0;x 1;x 2]};
    batesResult:@[batesFn;(trade;marketData;batesConfig);{x}];
    if[10h=type batesResult; :.limitcheck.errorRow[checkName;`bates;`heston;batesResult]];
    .limitcheck.resultRow[checkName;`bates;`heston;batesResult`unitPrice;hestonResult`unitPrice;absTol;relTol]
 };

/ --- Bates -> Merton ---

.modelcheck.checkBatesMertonLimit:{[trade;marketData;configDict]
    checkName:"BatesToMerton";
    mcConfig:.modelcheck.__getOr[configDict;`mcConfig;.modelcheck.__defaultMcConfig[]];
    absTol:.modelcheck.__getOr[configDict;`absoluteTolerance;1.0];
    relTol:.modelcheck.__getOr[configDict;`relativeTolerance;0.1];
    sigmaSquared:marketData[`volatility]*marketData`volatility;
    mertonParams:`volatility`jumpIntensity`jumpMean`jumpVolatility`riskFreeRate`dividendYield!(
        marketData`volatility;0.5;-0.1;0.3;marketData`riskFreeRate;marketData`dividendYield);
    mertonPrice:.merton.priceEuropeanSeries[trade`optionType;marketData`spot;trade`strike;trade`expiry;mertonParams;30];
    batesParams:`initialVariance`longRunVariance`meanReversion`volOfVol`correlation`jumpIntensity`jumpMean`jumpVolatility`riskFreeRate`dividendYield!(
        sigmaSquared;sigmaSquared;2.0;0.0;0.0;0.5;-0.1;0.3;marketData`riskFreeRate;marketData`dividendYield);
    batesConfig:`mcConfig`batesParams!(mcConfig;batesParams);
    batesFn:{.bates.priceEuropean[x 0;x 1;x 2]};
    batesResult:@[batesFn;(trade;marketData;batesConfig);{x}];
    if[10h=type batesResult; :.limitcheck.errorRow[checkName;`bates;`merton;batesResult]];
    .limitcheck.resultRow[checkName;`bates;`merton;batesResult`unitPrice;mertonPrice;absTol;relTol]
 };

/ --- Bates -> BS ---

.modelcheck.checkBatesBlackScholesLimit:{[trade;marketData;configDict]
    checkName:"BatesToBS";
    mcConfig:.modelcheck.__getOr[configDict;`mcConfig;.modelcheck.__defaultMcConfig[]];
    absTol:.modelcheck.__getOr[configDict;`absoluteTolerance;0.5];
    relTol:.modelcheck.__getOr[configDict;`relativeTolerance;0.05];
    bsPrice:.validation.blackScholesClosedForm[trade`optionType;marketData`spot;trade`strike;trade`expiry;marketData`riskFreeRate;marketData`dividendYield;marketData`volatility];
    sigmaSquared:marketData[`volatility]*marketData`volatility;
    batesParams:`initialVariance`longRunVariance`meanReversion`volOfVol`correlation`jumpIntensity`jumpMean`jumpVolatility`riskFreeRate`dividendYield!(
        sigmaSquared;sigmaSquared;2.0;0.0;0.0;0.0;0.0;0.0;marketData`riskFreeRate;marketData`dividendYield);
    batesConfig:`mcConfig`batesParams!(mcConfig;batesParams);
    batesFn:{.bates.priceEuropean[x 0;x 1;x 2]};
    batesResult:@[batesFn;(trade;marketData;batesConfig);{x}];
    if[10h=type batesResult; :.limitcheck.errorRow[checkName;`bates;`blackScholes;batesResult]];
    .limitcheck.resultRow[checkName;`bates;`blackScholes;batesResult`unitPrice;bsPrice;absTol;relTol]
 };

/ --- Local vol flat -> FDM ---

.modelcheck.checkLocalVolFlatLimit:{[trade;marketData;configDict]
    checkName:"LocalVolFlat";
    absTol:.modelcheck.__getOr[configDict;`absoluteTolerance;0.01];
    relTol:.modelcheck.__getOr[configDict;`relativeTolerance;0.001];
    bsModel:.model.createBlackScholesModel[];
    fdmConfig:.config.defaultPricingConfig[];
    constVol:marketData`volatility;
    flatVolFn:{[cv;spotValue;timePoint] cv}[constVol;;];
    lvModel:.model.createLocalVolatilityModel[];
    lvMkt:.market.createLocalVolatilityMarketData[marketData`underlying;marketData`spot;marketData`riskFreeRate;marketData`dividendYield;flatVolFn];
    engineFn:{.engine.priceOption[x 0;x 1;x 2;x 3]};
    flatResult:@[engineFn;(trade;marketData;bsModel;fdmConfig);{x}];
    if[10h=type flatResult; :.limitcheck.errorRow[checkName;`localVol;`fdm;flatResult]];
    lvResult:@[engineFn;(trade;lvMkt;lvModel;fdmConfig);{x}];
    if[10h=type lvResult; :.limitcheck.errorRow[checkName;`localVol;`fdm;lvResult]];
    .limitcheck.resultRow[checkName;`localVol;`fdm;lvResult`unitPrice;flatResult`unitPrice;absTol;relTol]
 };

/ --- SABR flat smile ---

.modelcheck.checkSabrFlatSmileLimit:{[smileInput;configDict]
    checkName:"SabrFlatSmile";
    absTol:.modelcheck.__getOr[configDict;`absoluteTolerance;0.01];
    relTol:.modelcheck.__getOr[configDict;`relativeTolerance;0.05];
    forwardValue:$[`forward in key smileInput;smileInput`forward;100f];
    expiryValue:$[`expiry in key smileInput;smileInput`expiry;1f];
    flatParams:`alpha`beta`rho`nu!(0.2;1.0;0.0;0.0001);
    atmVol:.sabr.impliedVolAtm[forwardValue;expiryValue;flatParams];
    otmVol:.sabr.impliedVol[forwardValue;120f;expiryValue;flatParams];
    .limitcheck.resultRow[checkName;`sabrOTM;`sabrATM;otmVol;atmVol;absTol;relTol]
 };

/ --- Barrier parity ---

.modelcheck.checkBarrierParity:{[trade;marketData;configDict]
    checkName:"BarrierParity";
    absTol:.modelcheck.__getOr[configDict;`absoluteTolerance;0.01];
    relTol:.modelcheck.__getOr[configDict;`relativeTolerance;0.001];
    bsModel:.model.createBlackScholesModel[];
    fdmConfig:.config.defaultPricingConfig[];
    / Vanilla
    vanillaTrade:trade;
    vanillaTrade[`barrierType]:`none;
    vanillaTrade[`barrierLevel]:0Nf;
    vanillaTrade[`rebate]:0f;
    engineFn:{.engine.priceOption[x 0;x 1;x 2;x 3]};
    vanillaResult:@[engineFn;(vanillaTrade;marketData;bsModel;fdmConfig);{x}];
    if[10h=type vanillaResult; :.limitcheck.errorRow[checkName;`barrierSum;`vanilla;vanillaResult]];
    / KO
    koTrade:trade;
    koTrade[`barrierType]:`upAndOut;
    koTrade[`barrierLevel]:120f;
    koTrade[`rebate]:0f;
    koResult:@[engineFn;(koTrade;marketData;bsModel;fdmConfig);{x}];
    if[10h=type koResult; :.limitcheck.errorRow[checkName;`barrierSum;`vanilla;koResult]];
    / KI
    kiTrade:trade;
    kiTrade[`barrierType]:`upAndIn;
    kiTrade[`barrierLevel]:120f;
    kiTrade[`rebate]:0f;
    kiResult:@[engineFn;(kiTrade;marketData;bsModel;fdmConfig);{x}];
    if[10h=type kiResult; :.limitcheck.errorRow[checkName;`barrierSum;`vanilla;kiResult]];
    barrierSum:koResult[`unitPrice]+kiResult`unitPrice;
    .limitcheck.resultRow[checkName;`barrierSum;`vanilla;barrierSum;vanillaResult`unitPrice;absTol;relTol]
 };

/ --- American-European consistency ---

.modelcheck.checkAmericanEuropeanConsistency:{[trade;marketData;configDict]
    checkName:"AmericanGeqEuropean";
    bsModel:.model.createBlackScholesModel[];
    fdmConfig:.config.defaultPricingConfig[];
    eurTrade:@[trade;`exerciseStyle;:;`european];
    engineFn:{.engine.priceOption[x 0;x 1;x 2;x 3]};
    eurResult:@[engineFn;(eurTrade;marketData;bsModel;fdmConfig);{x}];
    if[10h=type eurResult; :.limitcheck.errorRow[checkName;`american;`european;eurResult]];
    amerTrade:@[trade;`exerciseStyle;:;`american];
    amerResult:@[engineFn;(amerTrade;marketData;bsModel;fdmConfig);{x}];
    if[10h=type amerResult; :.limitcheck.errorRow[checkName;`american;`european;amerResult]];
    passedVal:amerResult[`unitPrice]>=eurResult[`unitPrice]-0.01;
    absDiff:amerResult[`unitPrice]-eurResult`unitPrice;
    `checkName`leftModel`rightModel`leftPrice`rightPrice`absoluteDifference`relativeDifference`absoluteTolerance`relativeTolerance`passed`status`errorMessage!(
        checkName;`american;`european;amerResult`unitPrice;eurResult`unitPrice;abs absDiff;0f;0f;0f;passedVal;`OK;"")
 };

/ --- Run all core checks ---

.modelcheck.runCoreModelLimitChecks:{[configDict]
    trade:.modelcheck.__defaultTrade[];
    mkt:.modelcheck.__defaultMarket[];
    resultTable:();
    / BS limits
    resultTable:resultTable,enlist .modelcheck.checkHestonBlackScholesLimit[trade;mkt;configDict];
    resultTable:resultTable,enlist .modelcheck.checkMertonBlackScholesLimit[trade;mkt;configDict];
    resultTable:resultTable,enlist .modelcheck.checkBatesBlackScholesLimit[trade;mkt;configDict];
    / Model hierarchy
    resultTable:resultTable,enlist .modelcheck.checkBatesHestonLimit[trade;mkt;configDict];
    resultTable:resultTable,enlist .modelcheck.checkBatesMertonLimit[trade;mkt;configDict];
    / FDM
    resultTable:resultTable,enlist .modelcheck.checkLocalVolFlatLimit[trade;mkt;configDict];
    / SABR
    resultTable:resultTable,enlist .modelcheck.checkSabrFlatSmileLimit[()!();configDict];
    / Barrier parity
    barrierTrade:trade,`barrierType`barrierLevel`rebate!(`upAndOut;120f;0f);
    resultTable:resultTable,enlist .modelcheck.checkBarrierParity[barrierTrade;mkt;configDict];
    / American >= European
    putTrade:@[trade;`optionType;:;`put];
    putTrade:putTrade,`barrierType`barrierLevel`rebate!(`none;0Nf;0f);
    divMkt:@[mkt;`dividendYield;:;0.03];
    resultTable:resultTable,enlist .modelcheck.checkAmericanEuropeanConsistency[putTrade;divMkt;configDict];
    resultTable
 };

/ --- Consistency report ---

.modelcheck.modelConsistencyReport:{[checkResultTable]
    resultRows:();
    loopIdx:0;
    while[loopIdx<count checkResultTable;
          rowData:checkResultTable loopIdx;
          severityVal:`pass;
          if[not rowData`passed;
             absDiff:rowData`absoluteDifference;
             absTol:rowData`absoluteTolerance;
             severityVal:$[absTol>0f;
                           $[absDiff<=2f*absTol;`warning;`fail];
                           `fail]];
          if[rowData[`status]~`ERROR; severityVal:`error];
          resultRows:resultRows,enlist `checkName`modelFamily`leftModel`rightModel`absoluteDifference`relativeDifference`passed`severity`status`errorMessage!(
              rowData`checkName;`limitCheck;rowData`leftModel;rowData`rightModel;
              rowData`absoluteDifference;rowData`relativeDifference;rowData`passed;severityVal;
              rowData`status;rowData`errorMessage);
          loopIdx+:1];
    resultRows
 };
