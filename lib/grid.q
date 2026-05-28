/ grid.q — spot and time grid construction

.grid.buildSpotGrid:{[sMin;sMax;n]
    sMin + ((sMax-sMin)%n) * til n+1
 };

.grid.buildTimeGrid:{[expiry;n]
    (expiry%n) * til n+1
 };

.grid.buildFiniteDifferenceGrid:{[trade;marketData;config]
    sMin:config`minimumSpot;
    sMax:config`maximumSpot;
    nSpot:config`numberOfSpotSteps;
    nTime:config`numberOfTimeSteps;
    spot:marketData`spot;
    strike:trade`strike;
    if[sMax <= spot; '"maximumSpot (",string[sMax],") must be greater than spot (",string[spot],")"];
    if[sMax <= strike; '"maximumSpot (",string[sMax],") must be greater than strike (",string[strike],")"];
    dS:(sMax-sMin) % nSpot;
    dt:trade[`expiry] % nTime;
    `spotGrid`timeGrid`spotStep`timeStep!(
        .grid.buildSpotGrid[sMin;sMax;nSpot];
        .grid.buildTimeGrid[trade`expiry;nTime];
        dS; dt)
 };
