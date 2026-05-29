/ commoditySpread.q - spread option pricing (Kirk approximation) (v0.32)

/ params dict keys: fwd1, fwd2, strike, expiry, vol1, vol2, correlation, riskFreeRate
.commodity.spread.validateInputs:{[optType;params]
    if[params[`fwd1]<=0f; '"Non-positive forward 1"];
    if[params[`fwd2]<=0f; '"Non-positive forward 2"];
    if[params[`expiry]<=0f; '"Non-positive expiry"];
    if[params[`vol1]<=0f; '"Non-positive vol1"];
    if[params[`vol2]<=0f; '"Non-positive vol2"];
    corrVal:params`correlation;
    if[(corrVal< -1f) or corrVal>1f; '"Correlation out of range"];
    if[not optType in `call`put; '"Invalid option type"];
 };

.commodity.spread.spreadPayoff:{[fwd1;fwd2;strikePrice;optType]
    spreadVal:fwd1-fwd2;
    $[optType=`call;0f|spreadVal-strikePrice;0f|strikePrice-spreadVal]
 };

/ Kirk approximation for spread option
/ params dict keys: fwd1, fwd2, strike, expiry, vol1, vol2, correlation, riskFreeRate
.commodity.spread.kirkPrice:{[optType;params]
    .commodity.spread.validateInputs[optType;params];
    fwd1:params`fwd1; fwd2:params`fwd2;
    strikePrice:params`strike; expiryVal:params`expiry;
    vol1:params`vol1; vol2:params`vol2;
    corrVal:params`correlation; rateVal:params`riskFreeRate;
    adjStrike:fwd2+strikePrice;
    if[adjStrike<=0f; '"Non-positive adjusted strike"];
    f2Ratio:fwd2%adjStrike;
    vol2Adj:vol2*f2Ratio;
    twoCorr:2f*corrVal*vol1*vol2Adj;
    kirkVol:sqrt (vol1*vol1)+(vol2Adj*vol2Adj)-twoCorr;
    .commodity.black76.price[optType;fwd1;adjStrike;expiryVal;kirkVol;rateVal]
 };

/ Calendar spread option (long near, short far)
.commodity.spread.calendarSpreadPrice:{[optType;params]
    .commodity.spread.kirkPrice[optType;params]
 };

/ -----------------------------------------------------------------------------
/ Margrabe (1978) exchange option (v0.49)
/ -----------------------------------------------------------------------------
/ Closed-form price of the option to exchange asset 2 for asset 1, i.e. payoff
/ max(F1 - F2, 0) for a `call and max(F2 - F1, 0) for a `put. This is the
/ zero-strike spread option and is exact (no Kirk-style strike adjustment).
/ params dict keys: fwd1, fwd2, expiry, vol1, vol2, correlation, riskFreeRate
/ (note: no strike key - the exchange option strike is identically zero).
.commodity.spread.__validateMargrabe:{[optType;params]
    if[params[`fwd1]<=0f; '"Non-positive forward 1"];
    if[params[`fwd2]<=0f; '"Non-positive forward 2"];
    if[params[`expiry]<=0f; '"Non-positive expiry"];
    if[params[`vol1]<=0f; '"Non-positive vol1"];
    if[params[`vol2]<=0f; '"Non-positive vol2"];
    corrVal:params`correlation;
    if[(corrVal< -1f) or corrVal>1f; '"Correlation out of range"];
    if[not optType in `call`put; '"Invalid option type"];
 };

.commodity.spread.margrabePrice:{[optType;params]
    .commodity.spread.__validateMargrabe[optType;params];
    fwd1:params`fwd1; fwd2:params`fwd2;
    expiryVal:params`expiry; vol1:params`vol1; vol2:params`vol2;
    corrVal:params`correlation; rateVal:params`riskFreeRate;
    exchVol:sqrt (vol1*vol1)+(vol2*vol2)-2f*corrVal*vol1*vol2;
    if[exchVol<=0f; '"Non-positive exchange volatility (perfectly correlated equal vols)"];
    sqrtT:sqrt expiryVal;
    d1Val:(log[fwd1%fwd2]+0.5*exchVol*exchVol*expiryVal)%(exchVol*sqrtT);
    d2Val:d1Val-exchVol*sqrtT;
    discFactor:exp neg rateVal*expiryVal;
    nd1:.validation.__normalCdf d1Val;
    nd2:.validation.__normalCdf d2Val;
    nNd1:.validation.__normalCdf neg d1Val;
    nNd2:.validation.__normalCdf neg d2Val;
    $[optType=`call;
        discFactor*(fwd1*nd1)-(fwd2*nd2);
        discFactor*(fwd2*nNd2)-(fwd1*nNd1)]
 };

/ -----------------------------------------------------------------------------
/ Monte Carlo spread option (v0.49)
/ -----------------------------------------------------------------------------
/ Two correlated lognormal forwards (martingales under their own measure, drift
/ -0.5 sigma^2 T), terminal payoff max(F1_T - F2_T - K, 0) for a `call. Uses the
/ shared MC normal generator + .montecarlo.priceFromPayoffs for discounting and
/ standard-error / confidence-interval reporting. Validates against Kirk and,
/ at strike 0, against the Margrabe closed form.
/ params dict keys: fwd1, fwd2, strike, expiry, vol1, vol2, correlation, riskFreeRate
/ mcConfig dict keys: pathCount, randomSeed; optional antithetic, confidenceLevel
.commodity.spread.spreadOptionMC:{[optType;params;mcConfig]
    .commodity.spread.validateInputs[optType;params];
    fwd1:params`fwd1; fwd2:params`fwd2; strikePrice:params`strike;
    expiryVal:params`expiry; vol1:params`vol1; vol2:params`vol2;
    corrVal:params`correlation; rateVal:params`riskFreeRate;
    pathCount:mcConfig`pathCount;
    seed:mcConfig`randomSeed;
    useAntithetic:$[`antithetic in key mcConfig; mcConfig`antithetic; 0b];
    confLevel:$[`confidenceLevel in key mcConfig; mcConfig`confidenceLevel; 0.95];
    if[pathCount<=0; '"spreadOptionMC pathCount must be positive"];
    nBase:$[useAntithetic; ceiling pathCount%2; pathCount];
    rawNormals:.montecarlo.__generateNormals[2*nBase;seed];
    eps1Base:nBase#rawNormals;
    eps2Base:nBase _ rawNormals;
    eps1:$[useAntithetic; eps1Base,neg eps1Base; eps1Base];
    eps2:$[useAntithetic; eps2Base,neg eps2Base; eps2Base];
    z1:eps1;
    z2:(corrVal*eps1)+(sqrt 1f-corrVal*corrVal)*eps2;
    sqrtT:sqrt expiryVal;
    f1Terminal:fwd1*exp (neg 0.5*vol1*vol1*expiryVal)+(vol1*sqrtT)*z1;
    f2Terminal:fwd2*exp (neg 0.5*vol2*vol2*expiryVal)+(vol2*sqrtT)*z2;
    spreadTerminal:f1Terminal-f2Terminal;
    payoffVector:$[optType=`call;0f|spreadTerminal-strikePrice;0f|strikePrice-spreadTerminal];
    .montecarlo.priceFromPayoffs[payoffVector;rateVal;expiryVal;confLevel]
 };
