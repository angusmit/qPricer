/ boundary.q — European option boundary conditions
/ Call: S=0 -> 0,  S=Smax -> Smax*exp(-q*tau) - K*exp(-r*tau)
/ Put:  S=0 -> K*exp(-r*tau),  S=Smax -> 0

.boundary.applyEuropeanBoundary:{[trade;marketData;grid;optVals;tau]
    strike:trade`strike;
    rDisc:exp neg marketData[`riskFreeRate]*tau;
    qDisc:exp neg marketData[`dividendYield]*tau;
    sMax:last grid`spotGrid;
    isCall:trade[`optionType]~`call;
    lo:$[isCall; 0f; strike*rDisc];
    hi:$[isCall; (sMax*qDisc)-strike*rDisc; 0f];
    @[@[optVals;0;:;lo];-1+count optVals;:;hi]
 };
