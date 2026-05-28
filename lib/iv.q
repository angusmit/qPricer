/ iv.q - Black-Scholes implied volatility solver (bisection)

.iv.__defaultSolverConfig:`lowerVolatilityBound`upperVolatilityBound`tolerance`maximumIterations!(0.0001;5.0;1e-8;100);

/ --- Public ---

.iv.impliedVolatility:{[optionType;marketPrice;spotPrice;strike;expiry;riskFreeRate;dividendYield]
    .iv.impliedVolatilityWithConfig[optionType;marketPrice;spotPrice;strike;expiry;riskFreeRate;dividendYield;.iv.__defaultSolverConfig]
 };

.iv.impliedVolatilityWithConfig:{[optionType;marketPrice;spotPrice;strike;expiry;riskFreeRate;dividendYield;solverConfig]
    .iv.__validateInputs[optionType;marketPrice;spotPrice;strike;expiry];
    .iv.__validateSolverConfig solverConfig;
    .iv.__checkNoArbitrageBounds[optionType;marketPrice;spotPrice;strike;expiry;riskFreeRate;dividendYield];
    lowerBound:solverConfig`lowerVolatilityBound;
    upperBound:solverConfig`upperVolatilityBound;
    tolerance:solverConfig`tolerance;
    maxIter:solverConfig`maximumIterations;
    priceAtLower:.validation.blackScholesClosedForm[optionType;spotPrice;strike;expiry;riskFreeRate;dividendYield;lowerBound];
    priceAtUpper:.validation.blackScholesClosedForm[optionType;spotPrice;strike;expiry;riskFreeRate;dividendYield;upperBound];
    if[marketPrice<priceAtLower; '"Implied volatility solver failed: market price below lower vol bound price"];
    if[marketPrice>priceAtUpper; '"Implied volatility solver failed: market price above upper vol bound price"];
    iterationCount:0;
    while[iterationCount<maxIter;
        candidateVolatility:0.5*lowerBound+upperBound;
        candidatePrice:.validation.blackScholesClosedForm[optionType;spotPrice;strike;expiry;riskFreeRate;dividendYield;candidateVolatility];
        priceDifference:candidatePrice-marketPrice;
        if[(abs priceDifference)<tolerance;
            :`impliedVolatility`iterations`pricingError`status`statusMessage!(
                candidateVolatility;iterationCount;priceDifference;`OK;"")];
        if[candidatePrice<marketPrice; lowerBound:candidateVolatility];
        if[candidatePrice>=marketPrice; upperBound:candidateVolatility];
        iterationCount+:1];
    finalVol:0.5*lowerBound+upperBound;
    finalPrice:.validation.blackScholesClosedForm[optionType;spotPrice;strike;expiry;riskFreeRate;dividendYield;finalVol];
    `impliedVolatility`iterations`pricingError`status`statusMessage!(
        finalVol;iterationCount;finalPrice-marketPrice;`OK;"converged within max iterations")
 };

.iv.calculateOptionChainImpliedVols:{[optionChain;spotPrice;riskFreeRate;dividendYield]
    numRows:count optionChain;
    resultList:();
    loopIdx:0;
    while[loopIdx<numRows;
        chainRow:optionChain loopIdx;
        ivResult:.[.iv.impliedVolatility;(chainRow`optionType;chainRow`marketPrice;spotPrice;chainRow`strike;chainRow`expiry;riskFreeRate;dividendYield);{x}];
        enrichedRow:chainRow;
        if[99h=type ivResult;
            enrichedRow:chainRow,`impliedVolatility`ivStatus`ivStatusMessage`ivPricingError!(
                ivResult`impliedVolatility;`OK;"";ivResult`pricingError)];
        if[10h=type ivResult;
            enrichedRow:chainRow,`impliedVolatility`ivStatus`ivStatusMessage`ivPricingError!(
                0Nf;`ERROR;ivResult;0Nf)];
        resultList:resultList,enlist enrichedRow;
        loopIdx+:1];
    resultList
 };

/ --- Internal ---

.iv.__validateInputs:{[optionType;marketPrice;spotPrice;strike;expiry]
    if[not optionType in `call`put; '"Implied vol: optionType must be `call or `put"];
    if[not marketPrice>0f; '"Implied vol: marketPrice must be positive"];
    if[not spotPrice>0f; '"Implied vol: spot must be positive"];
    if[not strike>0f; '"Implied vol: strike must be positive"];
    if[not expiry>0f; '"Implied vol: expiry must be positive"];
 };

.iv.__validateSolverConfig:{[solverConfig]
    if[not solverConfig[`lowerVolatilityBound]>0f; '"Implied vol: lowerVolatilityBound must be positive"];
    if[not solverConfig[`upperVolatilityBound]>solverConfig`lowerVolatilityBound;
        '"Implied vol: upperVolatilityBound must exceed lowerVolatilityBound"];
    if[not solverConfig[`maximumIterations]>0; '"Implied vol: maximumIterations must be positive"];
    if[not solverConfig[`tolerance]>0f; '"Implied vol: tolerance must be positive"];
 };

.iv.__checkNoArbitrageBounds:{[optionType;marketPrice;spotPrice;strike;expiry;riskFreeRate;dividendYield]
    discountFactor:exp neg riskFreeRate*expiry;
    dividendDiscount:exp neg dividendYield*expiry;
    priceBoundTolerance:1e-10;
    if[optionType~`call;
        lowerPriceBound:0f|(spotPrice*dividendDiscount)-strike*discountFactor;
        upperPriceBound:spotPrice*dividendDiscount;
        if[marketPrice<lowerPriceBound-priceBoundTolerance;
            '"Implied vol: call price ",string[marketPrice]," below no-arbitrage lower bound ",string lowerPriceBound];
        if[marketPrice>upperPriceBound+priceBoundTolerance;
            '"Implied vol: call price ",string[marketPrice]," above no-arbitrage upper bound ",string upperPriceBound]];
    if[optionType~`put;
        lowerPriceBound:0f|(strike*discountFactor)-spotPrice*dividendDiscount;
        upperPriceBound:strike*discountFactor;
        if[marketPrice<lowerPriceBound-priceBoundTolerance;
            '"Implied vol: put price ",string[marketPrice]," below no-arbitrage lower bound ",string lowerPriceBound];
        if[marketPrice>upperPriceBound+priceBoundTolerance;
            '"Implied vol: put price ",string[marketPrice]," above no-arbitrage upper bound ",string upperPriceBound]];
 };
