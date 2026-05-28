/ commodityBlack76.q - Black-76 commodity option pricing (v0.32)

.commodity.black76.validateInputs:{[fwdPrice;strikePrice;expiryVal;volVal;rateVal;optType]
    if[fwdPrice<=0f; '"Non-positive forward price"];
    if[strikePrice<=0f; '"Non-positive strike"];
    if[expiryVal<=0f; '"Non-positive expiry"];
    if[volVal<=0f; '"Non-positive volatility"];
    if[not optType in `call`put; '"Invalid option type: ",string optType];
 };

.commodity.black76.price:{[optType;fwdPrice;strikePrice;expiryVal;volVal;rateVal]
    .commodity.black76.validateInputs[fwdPrice;strikePrice;expiryVal;volVal;rateVal;optType];
    sqrtT:sqrt expiryVal;
    halfVolSq:0.5*volVal*volVal;
    d1Val:(log[fwdPrice%strikePrice]+halfVolSq*expiryVal)%(volVal*sqrtT);
    d2Val:d1Val-volVal*sqrtT;
    discFactor:exp neg rateVal*expiryVal;
    nd1:.validation.__normalCdf d1Val;
    nd2:.validation.__normalCdf d2Val;
    nNd1:.validation.__normalCdf neg d1Val;
    nNd2:.validation.__normalCdf neg d2Val;
    $[optType=`call;
        discFactor*(fwdPrice*nd1)-(strikePrice*nd2);
        discFactor*(strikePrice*nNd2)-(fwdPrice*nNd1)]
 };

.commodity.black76.greeks:{[optType;fwdPrice;strikePrice;expiryVal;volVal;rateVal]
    sqrtT:sqrt expiryVal;
    halfVolSq:0.5*volVal*volVal;
    d1Val:(log[fwdPrice%strikePrice]+halfVolSq*expiryVal)%(volVal*sqrtT);
    d2Val:d1Val-volVal*sqrtT;
    discFactor:exp neg rateVal*expiryVal;
    nd1:.validation.__normalCdf d1Val;
    nNd1:.validation.__normalCdf neg d1Val;
    pdfD1:exp[neg 0.5*d1Val*d1Val]%sqrt 2f*acos neg 1f;
    deltaVal:$[optType=`call;discFactor*nd1;neg discFactor*nNd1];
    gammaVal:discFactor*pdfD1%(fwdPrice*volVal*sqrtT);
    vegaVal:fwdPrice*discFactor*pdfD1*sqrtT%100f;
    thetaVal:neg fwdPrice*discFactor*pdfD1*volVal%(2f*sqrtT*365f);
    optPx:.commodity.black76.price[optType;fwdPrice;strikePrice;expiryVal;volVal;rateVal];
    rhoVal:neg expiryVal*discFactor*optPx%100f;
    `delta`gamma`vega`theta`rho!(deltaVal;gammaVal;vegaVal;thetaVal;rhoVal)
 };

/ Implied vol via bisection
/ solverConfig dict must have keys: loVol, hiVol, tol
.commodity.black76.impliedVol:{[optType;mktPrice;fwdPrice;strikePrice;expiryVal;rateVal;solverConfig]
    if[mktPrice<=0f; '"Non-positive market price for IV"];
    loVol:solverConfig`loVol;
    hiVol:solverConfig`hiVol;
    tol:solverConfig`tol;
    iterCount:0;
    maxIter:100;
    midVol:0Nf;
    modelPx:0Nf;
    while[iterCount<maxIter;
        midVol:(loVol+hiVol)%2f;
        modelPx:.commodity.black76.price[optType;fwdPrice;strikePrice;expiryVal;midVol;rateVal];
        diff:modelPx-mktPrice;
        if[(abs diff)<tol;
            :`impliedVol`iterations`modelPrice`marketPrice`status`errorMessage!(midVol;iterCount;modelPx;mktPrice;`OK;"")];
        if[diff>0f; hiVol:midVol];
        if[diff<0f; loVol:midVol];
        iterCount+:1];
    `impliedVol`iterations`modelPrice`marketPrice`status`errorMessage!(midVol;iterCount;modelPx;mktPrice;`warning;"max iterations reached")
 };
