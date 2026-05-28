/ schwartz2.q - Schwartz two-factor commodity log mean-reversion model (v0.34)
/ Risk-neutral dynamics:
/     d X_t = -kappa * X_t dt + sigmaX * dW1_t          (short-term mean-reverting)
/     d Y_t = muY dt + sigmaY * dW2_t                    (long-term factor)
/     corr(dW1, dW2) = rho
/     log S_t = X_t + Y_t
/ v0.34 treats all params as risk-neutral; market price of risk is deferred
/ to calibration. Initial factors shortFactor0 and longFactor0 are function
/ arguments, not part of the params dict, matching v0.33 conventions.
/ Public:  validateParams, factorMoments, logPriceMoments, futuresPrice,
/          futuresCurve, simulateFactors, simulatePaths, pathDiagnostics,
/          europeanOptionPrice, europeanOptionPriceMC.
/ Private: __shortFactorVariance, __crossCovariance, __exactShortStep,
/          __longStep, __correlatedPairs.

/ Threshold below which |kappa * tau| triggers Taylor branches.
.commodity.schwartz2.__smallKappaTau:1e-6;

.commodity.schwartz2.validateParams:{[params]
    requiredKeys:`meanReversionSpeed`shortVolatility`longVolatility`longDrift`correlation;
    paramKeys:key params;
    missingKeys:requiredKeys where not requiredKeys in paramKeys;
    if[0<count missingKeys;
        missingStr:", " sv string missingKeys;
        '"Schwartz2 params missing: ",missingStr];
    if[params[`meanReversionSpeed]<=0f; '"Schwartz2 meanReversionSpeed must be positive"];
    if[params[`shortVolatility]<0f; '"Schwartz2 shortVolatility must be non-negative"];
    if[params[`longVolatility]<0f; '"Schwartz2 longVolatility must be non-negative"];
    corrVal:params`correlation;
    if[corrVal<=-1f; '"Schwartz2 correlation must be greater than -1"];
    if[corrVal>=1f; '"Schwartz2 correlation must be less than 1"];
 };

/ Var[X_tau] = sigmaX^2 * (1 - exp(-2 kappa tau)) / (2 kappa).
/ Taylor branch for small |kappa*tau|:
/   sigmaX^2 * tau * (1 - kappa*tau + (2/3)*(kappa*tau)^2).
.commodity.schwartz2.__shortFactorVariance:{[params;tau]
    kappaVal:params`meanReversionSpeed;
    sigmaXVal:params`shortVolatility;
    sigmaXSq:sigmaXVal*sigmaXVal;
    kappaTau:kappaVal*tau;
    if[(abs kappaTau)<.commodity.schwartz2.__smallKappaTau;
        taylorFactor:(1f-kappaTau)+(2f%3f)*kappaTau*kappaTau;
        :sigmaXSq*tau*taylorFactor];
    sigmaXSq*(1f-exp neg 2f*kappaTau)%(2f*kappaVal)
 };

/ Cov[X_tau, Y_tau] = rho * sigmaX * sigmaY * (1 - exp(-kappa tau)) / kappa.
/ Taylor branch for small |kappa*tau|:
/   rho * sigmaX * sigmaY * tau * (1 - (kappa*tau)/2 + (kappa*tau)^2/6).
.commodity.schwartz2.__crossCovariance:{[params;tau]
    kappaVal:params`meanReversionSpeed;
    sigmaXVal:params`shortVolatility;
    sigmaYVal:params`longVolatility;
    rhoVal:params`correlation;
    kappaTau:kappaVal*tau;
    leadingTerm:rhoVal*sigmaXVal*sigmaYVal;
    if[(abs kappaTau)<.commodity.schwartz2.__smallKappaTau;
        taylorFactor:(1f-0.5*kappaTau)+(kappaTau*kappaTau)%6f;
        :leadingTerm*tau*taylorFactor];
    leadingTerm*(1f-exp neg kappaTau)%kappaVal
 };

/ Factor moments: per-factor means and variances plus cross-covariance.
.commodity.schwartz2.factorMoments:{[shortFactor0;longFactor0;params;tau]
    .commodity.schwartz2.validateParams params;
    if[tau<0f; '"Schwartz2 tau must be non-negative"];
    kappaVal:params`meanReversionSpeed;
    sigmaYVal:params`longVolatility;
    muYVal:params`longDrift;
    decayVal:exp neg kappaVal*tau;
    shortMeanVal:shortFactor0*decayVal;
    longMeanVal:longFactor0+muYVal*tau;
    shortVarianceVal:.commodity.schwartz2.__shortFactorVariance[params;tau];
    longVarianceVal:sigmaYVal*sigmaYVal*tau;
    crossCovarianceVal:.commodity.schwartz2.__crossCovariance[params;tau];
    `shortMean`longMean`shortVariance`longVariance`crossCovariance!(
        shortMeanVal;longMeanVal;shortVarianceVal;longVarianceVal;crossCovarianceVal)
 };

/ Combined log-price distribution: log S_tau ~ Normal(m, v).
.commodity.schwartz2.logPriceMoments:{[shortFactor0;longFactor0;params;tau]
    fm:.commodity.schwartz2.factorMoments[shortFactor0;longFactor0;params;tau];
    meanVal:fm[`shortMean]+fm`longMean;
    varianceVal:fm[`shortVariance]+fm[`longVariance]+2f*fm`crossCovariance;
    `mean`variance!(meanVal;varianceVal)
 };

/ Forward / futures price F(0,tau) = E[exp(X_tau + Y_tau)] = exp(m + v/2).
.commodity.schwartz2.futuresPrice:{[shortFactor0;longFactor0;params;tau]
    lpm:.commodity.schwartz2.logPriceMoments[shortFactor0;longFactor0;params;tau];
    exp lpm[`mean]+0.5*lpm`variance
 };

.commodity.schwartz2.futuresCurve:{[shortFactor0;longFactor0;params;tenorVector]
    .commodity.schwartz2.validateParams params;
    .commodity.schwartz2.futuresPrice[shortFactor0;longFactor0;params;] each tenorVector
 };

/ Exact OU step for the short factor.
.commodity.schwartz2.__exactShortStep:{[currentShortFactor;params;dtVal;normalDraws]
    kappaVal:params`meanReversionSpeed;
    decayVal:exp neg kappaVal*dtVal;
    stepVarianceVal:.commodity.schwartz2.__shortFactorVariance[params;dtVal];
    stepStdVal:sqrt stepVarianceVal;
    (decayVal*currentShortFactor)+stepStdVal*normalDraws
 };

/ Arithmetic Brownian step for the long factor.
.commodity.schwartz2.__longStep:{[currentLongFactor;params;dtVal;normalDraws]
    sigmaYVal:params`longVolatility;
    muYVal:params`longDrift;
    currentLongFactor+(muYVal*dtVal)+sigmaYVal*(sqrt dtVal)*normalDraws
 };

/ Correlated normal pair generator with antithetic support.
/ Returns dict `z1`z2 of pathCount x stepCount matrices with corr(Z1,Z2) = rho.
/ Antithetic preserved by negating BOTH Z1 and Z2 on the back half of paths.
.commodity.schwartz2.__correlatedPairs:{[pathCount;stepCount;correlation;randomSeed;antithetic]
    if[not antithetic;
        :.heston.generateCorrelatedBrownian[pathCount;stepCount;correlation;randomSeed]];
    halfPaths:ceiling pathCount%2;
    basePair:.heston.generateCorrelatedBrownian[halfPaths;stepCount;correlation;randomSeed];
    z1Base:basePair`z1;
    z2Base:basePair`z2;
    z1Full:pathCount#(z1Base,neg z1Base);
    z2Full:pathCount#(z2Base,neg z2Base);
    `z1`z2!(z1Full;z2Full)
 };

/ Simulate factor paths. Returns dict `shortPaths`longPaths, each pathCount x stepCount
/ (excluding t=0), matching the .montecarlo.simulateGBMPaths convention.
.commodity.schwartz2.simulateFactors:{[shortFactor0;longFactor0;params;expiry;mcConfig]
    .commodity.schwartz2.validateParams params;
    .montecarlo.validateMcConfig mcConfig;
    if[expiry<=0f; '"Schwartz2 expiry must be positive"];
    pathCountValue:mcConfig`pathCount;
    timeStepCountValue:mcConfig`timeStepCount;
    randomSeedValue:mcConfig`randomSeed;
    useAntithetic:mcConfig`antithetic;
    rhoVal:params`correlation;
    dtVal:expiry%timeStepCountValue;
    brownianPair:.commodity.schwartz2.__correlatedPairs[pathCountValue;timeStepCountValue;rhoVal;randomSeedValue;useAntithetic];
    z1Matrix:brownianPair`z1;
    z2Matrix:brownianPair`z2;
    currentShortFactors:pathCountValue#shortFactor0;
    currentLongFactors:pathCountValue#longFactor0;
    shortRows:enlist currentShortFactors;
    longRows:enlist currentLongFactors;
    stepIdx:0;
    while[stepIdx<timeStepCountValue;
        z1AtStep:z1Matrix[;stepIdx];
        z2AtStep:z2Matrix[;stepIdx];
        currentShortFactors:.commodity.schwartz2.__exactShortStep[currentShortFactors;params;dtVal;z1AtStep];
        currentLongFactors:.commodity.schwartz2.__longStep[currentLongFactors;params;dtVal;z2AtStep];
        shortRows:shortRows,enlist currentShortFactors;
        longRows:longRows,enlist currentLongFactors;
        stepIdx+:1];
    `shortPaths`longPaths!(flip 1_shortRows;flip 1_longRows)
 };

/ Simulate spot price paths S_t = exp(X_t + Y_t).
.commodity.schwartz2.simulatePaths:{[shortFactor0;longFactor0;params;expiry;mcConfig]
    factors:.commodity.schwartz2.simulateFactors[shortFactor0;longFactor0;params;expiry;mcConfig];
    shortPathsMat:factors`shortPaths;
    longPathsMat:factors`longPaths;
    {[shortRow;longRow] exp shortRow+longRow}'[shortPathsMat;longPathsMat]
 };

/ Path diagnostics wrapper using .pathdiag.summary on the spot path matrix.
.commodity.schwartz2.pathDiagnostics:{[pathMatrix]
    .pathdiag.summary[pathMatrix;1f]
 };

/ Closed-form European option: S_T is exactly lognormal under this Gaussian
/ two-factor setup, so reuse Black-76 with the implied forward and effective vol.
.commodity.schwartz2.europeanOptionPrice:{[optType;shortFactor0;longFactor0;strikePrice;expiry;riskFreeRate;params]
    .commodity.schwartz2.validateParams params;
    if[not optType in `call`put; '"Schwartz2 option type must be call or put"];
    if[strikePrice<=0f; '"Schwartz2 strike must be positive"];
    if[expiry<=0f; '"Schwartz2 expiry must be positive"];
    lpm:.commodity.schwartz2.logPriceMoments[shortFactor0;longFactor0;params;expiry];
    meanVal:lpm`mean;
    varianceVal:lpm`variance;
    if[varianceVal<0f; '"Schwartz2 log variance went negative; check parameters"];
    fwdPrice:exp meanVal+0.5*varianceVal;
    effectiveVol:sqrt varianceVal%expiry;
    .commodity.black76.price[optType;fwdPrice;strikePrice;expiry;effectiveVol;riskFreeRate]
 };

/ Monte Carlo cross-check.
.commodity.schwartz2.europeanOptionPriceMC:{[optType;shortFactor0;longFactor0;strikePrice;expiry;riskFreeRate;params;mcConfig]
    if[not optType in `call`put; '"Schwartz2 option type must be call or put"];
    if[strikePrice<=0f; '"Schwartz2 strike must be positive"];
    pathMatrix:.commodity.schwartz2.simulatePaths[shortFactor0;longFactor0;params;expiry;mcConfig];
    terminalSpots:last each pathMatrix;
    payoffVector:$[optType=`call;0f|terminalSpots-strikePrice;0f|strikePrice-terminalSpots];
    .montecarlo.priceFromPayoffs[payoffVector;riskFreeRate;expiry;mcConfig`confidenceLevel]
 };
