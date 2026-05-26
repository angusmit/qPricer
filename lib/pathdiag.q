/ pathdiag.q - path diagnostic utilities for Monte Carlo simulations (v0.18)

.pathdiag.pathMinimum:{[pathMatrix] min each pathMatrix};
.pathdiag.pathMaximum:{[pathMatrix] max each pathMatrix};
.pathdiag.pathFinal:{[pathMatrix] last each pathMatrix};
.pathdiag.pathAverage:{[pathMatrix] avg each pathMatrix};

.pathdiag.pathReturn:{[pathMatrix]
    finalSpots:last each pathMatrix;
    initialSpots:first each pathMatrix;
    finalSpots%initialSpots
 };

.pathdiag.pathRealisedVol:{[pathMatrix;expiry]
    logReturns:{1_ deltas log x} each pathMatrix;
    stepVols:dev each logReturns;
    stepsPerPath:count pathMatrix 0;
    dtVal:expiry%stepsPerPath;
    stepVols%sqrt dtVal
 };

.pathdiag.summary:{[pathMatrix;expiry]
    finalVector:.pathdiag.pathFinal pathMatrix;
    minVector:.pathdiag.pathMinimum pathMatrix;
    maxVector:.pathdiag.pathMaximum pathMatrix;
    returnVector:.pathdiag.pathReturn pathMatrix;
    realisedVolVector:.pathdiag.pathRealisedVol[pathMatrix;expiry];
    `pathCount`timeStepCount`averageFinal`averageMinimum`averageMaximum`averageReturn`averageRealisedVol!(
        count pathMatrix;
        count pathMatrix 0;
        avg finalVector;
        avg minVector;
        avg maxVector;
        avg returnVector;
        avg realisedVolVector)
 };
