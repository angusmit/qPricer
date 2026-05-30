/ variance.q - variance reduction utilities for Monte Carlo (v0.19)

/ --- Antithetic ---

.variance.applyAntithetic:{[normalMatrix]
    halfPaths:ceiling (count normalMatrix)%2;
    basePaths:halfPaths#normalMatrix;
    antiPaths:neg each basePaths;
    (count normalMatrix)#basePaths,antiPaths
 };

/ --- Moment matching ---

.variance.applyMomentMatching:{[normalMatrix]
    stepsPerPath:count normalMatrix 0;
    allNormals:raze normalMatrix;
    meanVal:avg allNormals;
    stdVal:dev allNormals;
    if[stdVal<=0f; :normalMatrix];
    adjustedFlat:(allNormals-meanVal)%stdVal;
    stepsPerPath cut adjustedFlat
 };

/ --- Variance/SE estimation ---

.variance.estimateVariance:{[payoffVector] var payoffVector};

.variance.standardError:{[payoffVector]
    (dev payoffVector)%sqrt `float$count payoffVector
 };

.variance.confidenceInterval:{[priceVal;standardError;confidenceLevel]
    zScore:.montecarlo.__zScoreForConfidence confidenceLevel;
    `lowerConfidence`upperConfidence!(priceVal-zScore*standardError;priceVal+zScore*standardError)
 };

/ --- Control variate ---

.variance.controlVariateBeta:{[targetPayoff;controlPayoff]
    covValue:cov[targetPayoff;controlPayoff];
    varValue:var controlPayoff;
    if[varValue<=0f; :0f];
    covValue%varValue
 };

.variance.controlVariateAdjust:{[targetPayoff;controlPayoff;controlExpectedPayoff]
    betaValue:.variance.controlVariateBeta[targetPayoff;controlPayoff];
    targetPayoff-betaValue*(controlPayoff-controlExpectedPayoff)
 };

/ --- Comparison ---

.variance.compareStandardErrors:{[plainSE;reducedSE]
    reductionRatio:$[plainSE>0f;reducedSE%plainSE;0Nf];
    `plainSE`reducedSE`reductionRatio`varianceReductionFactor!(
        plainSE;reducedSE;reductionRatio;$[reductionRatio>0f;1f%reductionRatio*reductionRatio;0Nf])
 };
