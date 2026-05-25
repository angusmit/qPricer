/ iv.q - Black-Scholes implied volatility solver (bisection)
/ Recovers the volatility that reproduces a given market price.

.iv.__defaultSolverConfig:`lowerVolatilityBound`upperVolatilityBound`tolerance`maximumIterations!(0.0001;5.0;1e-8;100);

/ --- Public ---

.iv.impliedVolatility:{[optionType;marketPrice;spot;strike;expiry;riskFreeRate;dividendYield]
    .iv.impliedVolatilityWithConfig[optionType;marketPrice;spot;strike;expiry;riskFreeRate;dividendYield;.iv.__defaultSolverConfig]
 };

.iv.impliedVolatilityWithConfig:{[optionType;marketPrice;spot;strike;expiry;riskFreeRate;dividendYield;solverConfig]
    .iv.__validateInputs[optionType;marketPrice;spot;strike;expiry];
    lowerBound:solverConfig`lowerVolatilityBound;
    upperBound:solverConfig`upperVolatilityBound;
    tolerance:solverConfig`tolerance;
    maxIter:solverConfig`maximumIterations;
    priceAtLower:.validation.blackScholesClosedForm[optionType;spot;strike;expiry;riskFreeRate;dividendYield;lowerBound];
    priceAtUpper:.validation.blackScholesClosedForm[optionType;spot;strike;expiry;riskFreeRate;dividendYield;upperBound];
    if[marketPrice<priceAtLower; '"Implied volatility solver failed: market price below lower vol bound price"];
    if[marketPrice>priceAtUpper; '"Implied volatility solver failed: market price above upper vol bound price"];
    / Bisection loop
    iterationCount:0;
    while[iterationCount<maxIter;
        candidateVolatility:0.5*lowerBound+upperBound;
        candidatePrice:.validation.blackScholesClosedForm[optionType;spot;strike;expiry;riskFreeRate;dividendYield;candidateVolatility];
        priceDifference:candidatePrice-marketPrice;
        if[(abs priceDifference)<tolerance;
            :`impliedVolatility`iterations`pricingError`status`statusMessage!(
                candidateVolatility;iterationCount;priceDifference;`OK;"")];
        if[candidatePrice<marketPrice; lowerBound:candidateVolatility];
        if[candidatePrice>=marketPrice; upperBound:candidateVolatility];
        iterationCount+:1];
    / Max iterations reached - return best estimate
    finalVol:0.5*lowerBound+upperBound;
    finalPrice:.validation.blackScholesClosedForm[optionType;spot;strike;expiry;riskFreeRate;dividendYield;finalVol];
    `impliedVolatility`iterations`pricingError`status`statusMessage!(
        finalVol;iterationCount;finalPrice-marketPrice;`OK;"converged within max iterations")
 };

.iv.calculateOptionChainImpliedVols:{[optionChain;spot;riskFreeRate;dividendYield]
    numRows:count optionChain;
    resultList:();
    loopIdx:0;
    while[loopIdx<numRows;
        chainRow:optionChain loopIdx;
        ivResult:.[.iv.impliedVolatility;(chainRow`optionType;chainRow`marketPrice;spot;chainRow`strike;chainRow`expiry;riskFreeRate;dividendYield);{x}];
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

.iv.__validateInputs:{[optionType;marketPrice;spot;strike;expiry]
    if[not optionType in `call`put; '"Implied vol: optionType must be `call or `put"];
    if[not marketPrice>0f; '"Implied vol: marketPrice must be positive"];
    if[not spot>0f; '"Implied vol: spot must be positive"];
    if[not strike>0f; '"Implied vol: strike must be positive"];
    if[not expiry>0f; '"Implied vol: expiry must be positive"];
 };
