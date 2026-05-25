/ boundary.q - European boundary conditions + barrier knock-out
/ Call: S=0 -> 0,  S=Smax -> Smax*exp(-q*tau) - K*exp(-r*tau)
/ Put:  S=0 -> K*exp(-r*tau),  S=Smax -> 0
/ Barrier: zero out option values at/beyond barrier level

.boundary.applyEuropeanBoundary:{[trade;marketData;gridDict;optionValues;remainingTime]
    strike:trade`strike;
    rDisc:exp neg marketData[`riskFreeRate]*remainingTime;
    qDisc:exp neg marketData[`dividendYield]*remainingTime;
    sMax:last gridDict`spotGrid;
    isCall:trade[`optionType]~`call;
    lowBoundary:0f;
    if[not isCall; lowBoundary:strike*rDisc];
    highBoundary:0f;
    if[isCall; highBoundary:(sMax*qDisc)-strike*rDisc];
    optionValues:@[optionValues;0;:;lowBoundary];
    @[optionValues;(-1)+count optionValues;:;highBoundary]
 };

/ Apply knock-out barrier condition: zero out values at/beyond barrier
.boundary.applyBarrierCondition:{[trade;spotGrid;optionValues]
    barrierType:.product.getBarrierType trade;
    if[barrierType~`none; :optionValues];
    barrierLevel:.product.getBarrierLevel trade;
    / Up-and-out: zero where spot >= barrier
    if[barrierType~`upAndOut;
        knockOutMask:spotGrid>=barrierLevel;
        :@[optionValues;where knockOutMask;:;0f]];
    / Down-and-out: zero where spot <= barrier
    if[barrierType~`downAndOut;
        knockOutMask:spotGrid<=barrierLevel;
        :@[optionValues;where knockOutMask;:;0f]];
    optionValues
 };
