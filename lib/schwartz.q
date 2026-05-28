/ schwartz.q - Schwartz (1997) one-factor commodity log mean-reversion model (v0.33)
/ Risk-neutral dynamics of the log price:
/     d X_t = kappa * (theta - X_t) dt + sigma dW_t
/ where X_t is the log spot/futures proxy and F_t = exp X_t.
/ v0.33 treats longRunLogMean as risk-neutral; market price of risk is deferred
/ to calibration/future model extension.
/ Public:  validateParams, transitionMoments, stationaryVariance, futuresPrice,
/          futuresCurve, simulatePaths, europeanOptionPrice, europeanOptionPriceMC.
/ Private: __exactStep, __effectiveVariance, __effectiveBlackScholesVol.

/ Threshold below which |kappa * tau| triggers the Taylor branch in effectiveVariance.
.commodity.schwartz.__smallKappaTau:1e-6;

.commodity.schwartz.validateParams:{[params]
    requiredKeys:`meanReversionSpeed`longRunLogMean`volatility;
    paramKeys:key params;
    missingKeys:requiredKeys where not requiredKeys in paramKeys;
    if[0<count missingKeys;
        missingStr:", " sv string missingKeys;
        '"Schwartz params missing: ",missingStr];
    if[params[`meanReversionSpeed]<=0f; '"Schwartz meanReversionSpeed must be positive"];
    if[params[`volatility]<=0f; '"Schwartz volatility must be positive"];
 };

/ Integrated variance Var[X_tau | X_0] = sigma^2/(2 kappa) * (1 - exp(-2 kappa tau)).
/ Small-kappa Taylor:  sigma^2 * tau * (1 - kappa*tau + (2/3)*(kappa*tau)^2).
.commodity.schwartz.__effectiveVariance:{[params;tau]
    kappaVal:params`meanReversionSpeed;
    sigmaVal:params`volatility;
    sigmaSq:sigmaVal*sigmaVal;
    kappaTau:kappaVal*tau;
    if[(abs kappaTau)<.commodity.schwartz.__smallKappaTau;
        taylorFactor:(1f-kappaTau)+(2f%3f)*kappaTau*kappaTau;
        :sigmaSq*tau*taylorFactor];
    sigmaSq*(1f-exp neg 2f*kappaTau)%(2f*kappaVal)
 };

.commodity.schwartz.__effectiveBlackScholesVol:{[params;tau]
    sqrt .commodity.schwartz.__effectiveVariance[params;tau]%tau
 };

.commodity.schwartz.stationaryVariance:{[params]
    .commodity.schwartz.validateParams params;
    sigmaVal:params`volatility;
    kappaVal:params`meanReversionSpeed;
    sigmaVal*sigmaVal%(2f*kappaVal)
 };

/ Conditional mean and variance of X_tau given X_0 = x0.
.commodity.schwartz.transitionMoments:{[x0;params;tau]
    .commodity.schwartz.validateParams params;
    if[tau<0f; '"Schwartz tau must be non-negative"];
    kappaVal:params`meanReversionSpeed;
    thetaVal:params`longRunLogMean;
    decayVal:exp neg kappaVal*tau;
    meanVal:(decayVal*x0)+(1f-decayVal)*thetaVal;
    varianceVal:.commodity.schwartz.__effectiveVariance[params;tau];
    `mean`variance!(meanVal;varianceVal)
 };

/ Risk-neutral forward/futures price F(0,T) = E[exp X_T | X_0] = exp(m_T + v_T/2).
.commodity.schwartz.futuresPrice:{[x0;params;tau]
    momentsDict:.commodity.schwartz.transitionMoments[x0;params;tau];
    exp momentsDict[`mean]+0.5*momentsDict`variance
 };

.commodity.schwartz.futuresCurve:{[x0;params;tenorVector]
    .commodity.schwartz.validateParams params;
    .commodity.schwartz.futuresPrice[x0;params;] each tenorVector
 };

/ One exact OU transition step (vectorised over paths via normalDraws).
.commodity.schwartz.__exactStep:{[currentLogPrice;params;dtVal;normalDraws]
    kappaVal:params`meanReversionSpeed;
    thetaVal:params`longRunLogMean;
    decayVal:exp neg kappaVal*dtVal;
    stepVarianceVal:.commodity.schwartz.__effectiveVariance[params;dtVal];
    stepStdVal:sqrt stepVarianceVal;
    (decayVal*currentLogPrice)+((1f-decayVal)*thetaVal)+stepStdVal*normalDraws
 };

/ Returns pathCount paths of length timeStepCount (excluding the t=0 state),
/ matching the .montecarlo.simulateGBMPaths convention. Values are prices = exp X.
.commodity.schwartz.simulatePaths:{[x0;params;expiry;mcConfig]
    .commodity.schwartz.validateParams params;
    .montecarlo.validateMcConfig mcConfig;
    if[expiry<=0f; '"Schwartz expiry must be positive"];
    pathCountValue:mcConfig`pathCount;
    timeStepCountValue:mcConfig`timeStepCount;
    randomSeedValue:mcConfig`randomSeed;
    useAntithetic:mcConfig`antithetic;
    dtVal:expiry%timeStepCountValue;
    normalMatrix:.montecarlo.generateNormalPaths[pathCountValue;timeStepCountValue;randomSeedValue];
    if[useAntithetic; normalMatrix:.montecarlo.generateAntitheticNormals[pathCountValue;timeStepCountValue;randomSeedValue]];
    currentLogPrices:pathCountValue#x0;
    logPathRows:enlist currentLogPrices;
    stepIdx:0;
    while[stepIdx<timeStepCountValue;
        normalsAtStep:normalMatrix[;stepIdx];
        currentLogPrices:.commodity.schwartz.__exactStep[currentLogPrices;params;dtVal;normalsAtStep];
        logPathRows:logPathRows,enlist currentLogPrices;
        stepIdx+:1];
    flip exp 1_logPathRows
 };

/ Closed-form European option on spot under Schwartz one-factor.
/ S_T is lognormal with log-mean m_T, log-variance v_T from transitionMoments.
/ With forward F = exp(m_T + v_T/2) and integrated stddev s = sqrt v_T:
/   call = exp(-rT) * (F*N(d1) - K*N(d2)),  d1 = (log(F/K) + s^2/2)/s,  d2 = d1 - s.
.commodity.schwartz.europeanOptionPrice:{[optType;x0;strikePrice;expiry;riskFreeRate;params]
    .commodity.schwartz.validateParams params;
    if[not optType in `call`put; '"Schwartz option type must be call or put"];
    if[strikePrice<=0f; '"Schwartz strike must be positive"];
    if[expiry<=0f; '"Schwartz expiry must be positive"];
    fwdPx:.commodity.schwartz.futuresPrice[x0;params;expiry];
    sVal:sqrt .commodity.schwartz.__effectiveVariance[params;expiry];
    d1Val:((log fwdPx%strikePrice)+0.5*sVal*sVal)%sVal;
    d2Val:d1Val-sVal;
    discFactor:exp neg riskFreeRate*expiry;
    nd1:.validation.__normalCdf d1Val;
    nd2:.validation.__normalCdf d2Val;
    nNd1:.validation.__normalCdf neg d1Val;
    nNd2:.validation.__normalCdf neg d2Val;
    $[optType=`call;
        discFactor*(fwdPx*nd1)-strikePrice*nd2;
        discFactor*(strikePrice*nNd2)-fwdPx*nNd1]
 };

/ Monte Carlo cross-check for the European option.
.commodity.schwartz.europeanOptionPriceMC:{[optType;x0;strikePrice;expiry;riskFreeRate;params;mcConfig]
    if[not optType in `call`put; '"Schwartz option type must be call or put"];
    if[strikePrice<=0f; '"Schwartz strike must be positive"];
    pathMatrix:.commodity.schwartz.simulatePaths[x0;params;expiry;mcConfig];
    terminalSpots:last each pathMatrix;
    payoffVector:$[optType=`call;0f|terminalSpots-strikePrice;0f|strikePrice-terminalSpots];
    .montecarlo.priceFromPayoffs[payoffVector;riskFreeRate;expiry;mcConfig`confidenceLevel]
 };
