/ sabr.q - SABR implied volatility model and calibration (v0.23)

/ --- Validation ---

.sabr.validateParams:{[sabrParams]
    if[not sabrParams[`alpha]>0f; '"alpha must be positive"];
    if[sabrParams[`beta]<0f; '"beta must be >= 0"];
    if[sabrParams[`beta]>1f; '"beta must be <= 1"];
    if[sabrParams[`rho]<= -1f; '"rho must be > -1"];
    if[sabrParams[`rho]>=1f; '"rho must be < 1"];
    if[sabrParams[`nu]<0f; '"nu must be non-negative"];
 };

/ --- ATM implied vol ---

.sabr.impliedVolAtm:{[forwardValue;expiryValue;sabrParams]
    alphaVal:sabrParams`alpha;
    betaVal:sabrParams`beta;
    rhoVal:sabrParams`rho;
    nuVal:sabrParams`nu;
    oneMinusBeta:1f-betaVal;
    logF:log forwardValue;
    fPow:exp oneMinusBeta*logF;
    fPow2:fPow*fPow;
    omBeta2:oneMinusBeta*oneMinusBeta;
    alpha2:alphaVal*alphaVal;
    nu2:nuVal*nuVal;
    rho2:rhoVal*rhoVal;
    corrTerm1:omBeta2*alpha2%(24f*fPow2);
    corrTerm2:rhoVal*betaVal*nuVal*alphaVal%(4f*fPow);
    threeRho2:3f*rho2;
    corrTerm3:((2f-threeRho2)%24f)*nu2;
    corrFactor:1f+(corrTerm1+corrTerm2+corrTerm3)*expiryValue;
    (alphaVal%fPow)*corrFactor
 };

/ --- Non-ATM implied vol (Hagan approximation) ---

.sabr.impliedVol:{[forwardValue;strikeValue;expiryValue;sabrParams]
    .sabr.validateParams sabrParams;
    if[not forwardValue>0f; '"forward must be positive"];
    if[not strikeValue>0f; '"strike must be positive"];
    if[not expiryValue>0f; '"expiry must be positive"];
    alphaVal:sabrParams`alpha;
    betaVal:sabrParams`beta;
    rhoVal:sabrParams`rho;
    nuVal:sabrParams`nu;
    logFK:log forwardValue%strikeValue;
    / Near-ATM: use stable ATM formula
    if[(abs logFK)<1e-7; :.sabr.impliedVolAtm[forwardValue;expiryValue;sabrParams]];
    oneMinusBeta:1f-betaVal;
    fkProduct:forwardValue*strikeValue;
    logFkProduct:log fkProduct;
    fkPowHalf:exp 0.5*oneMinusBeta*logFkProduct;
    / z = (nu/alpha) * (FK)^((1-β)/2) * log(F/K)
    zVal:$[alphaVal>1e-15;(nuVal%alphaVal)*fkPowHalf*logFK;0f];
    / χ(z)
    twoRhoZ:2f*rhoVal*zVal;
    zSquared:zVal*zVal;
    sqrtArg:(1f-twoRhoZ)+zSquared;
    sqrtArg:0.0001|sqrtArg;
    sqrtTerm:sqrt sqrtArg;
    numerChi:(sqrtTerm+zVal)-rhoVal;
    denomChi:1f-rhoVal;
    chiVal:log numerChi%denomChi;
    zOverChi:$[(abs chiVal)>1e-12;zVal%chiVal;1f];
    / Denominator
    logFK2:logFK*logFK;
    logFK4:logFK2*logFK2;
    omBeta2:oneMinusBeta*oneMinusBeta;
    denomA:omBeta2%24f;
    denomB:omBeta2*omBeta2%1920f;
    denominator:fkPowHalf*(1f+(denomA*logFK2)+denomB*logFK4);
    / Correction factor
    fkPow:exp oneMinusBeta*logFkProduct;
    alpha2:alphaVal*alphaVal;
    nu2:nuVal*nuVal;
    rho2:rhoVal*rhoVal;
    corrA:omBeta2*alpha2%(24f*fkPow);
    corrB:rhoVal*betaVal*nuVal*alphaVal%(4f*fkPowHalf);
    threeRho2:3f*rho2;
    corrC:((2f-threeRho2)%24f)*nu2;
    corrFactor:1f+(corrA+corrB+corrC)*expiryValue;
    impliedVolVal:(alphaVal%denominator)*zOverChi*corrFactor;
    0.0001|impliedVolVal
 };

/ --- Smile and surface generation ---

.sabr.smileTable:{[forwardValue;strikeList;expiryValue;sabrParams]
    resultRows:();
    strikeIdx:0;
    while[strikeIdx<count strikeList;
        strikeVal:strikeList strikeIdx;
        ivResult:@[.sabr.impliedVol[forwardValue;strikeVal;expiryValue;];sabrParams;{x}];
        ivVal:$[10h=type ivResult;0Nf;ivResult];
        statusVal:$[10h=type ivResult;`ERROR;`OK];
        errMsg:$[10h=type ivResult;ivResult;""];
        resultRows:resultRows,enlist `forward`strike`expiry`impliedVolatility`alpha`beta`rho`nu`status`errorMessage!(
            forwardValue;strikeVal;expiryValue;ivVal;
            sabrParams`alpha;sabrParams`beta;sabrParams`rho;sabrParams`nu;
            statusVal;errMsg);
        strikeIdx+:1];
    resultRows
 };

.sabr.surfaceTable:{[forwardValue;strikeList;expiryList;sabrParams]
    resultRows:();
    expiryIdx:0;
    while[expiryIdx<count expiryList;
        expiryVal:expiryList expiryIdx;
        smileRows:.sabr.smileTable[forwardValue;strikeList;expiryVal;sabrParams];
        resultRows:resultRows,smileRows;
        expiryIdx+:1];
    resultRows
 };

/ --- Calibration ---

.sabr.validateMarketSmile:{[smileTable]
    if[0=count smileTable; '"Smile table is empty"];
    firstRow:smileTable 0;
    tableCols:key firstRow;
    if[not `marketImpliedVolatility in tableCols; '"Missing marketImpliedVolatility column"];
    if[not `strike in tableCols; '"Missing strike column"];
    if[not `expiry in tableCols; '"Missing expiry column"];
 };

.sabr.calibrationGrid:{[paramGridDict]
    alphaList:paramGridDict`alphaList;
    betaList:paramGridDict`betaList;
    rhoList:paramGridDict`rhoList;
    nuList:paramGridDict`nuList;
    gridRows:();
    aIdx:0;
    while[aIdx<count alphaList;
        bIdx:0;
        while[bIdx<count betaList;
            rIdx:0;
            while[rIdx<count rhoList;
                nIdx:0;
                while[nIdx<count nuList;
                    gridRows:gridRows,enlist `alpha`beta`rho`nu!(alphaList aIdx;betaList bIdx;rhoList rIdx;nuList nIdx);
                    nIdx+:1];
                rIdx+:1];
            bIdx+:1];
        aIdx+:1];
    gridRows
 };

.sabr.calibrateSmile:{[smileTable;forwardValue;sabrGrid;configDict]
    .sabr.validateMarketSmile smileTable;
    calibrationTable:();
    gridIdx:0;
    while[gridIdx<count sabrGrid;
        paramRow:sabrGrid gridIdx;
        sabrParams:`alpha`beta`rho`nu!(paramRow`alpha;paramRow`beta;paramRow`rho;paramRow`nu);
        modelVols:();
        marketVols:();
        weights:();
        failedCount:0;
        smileIdx:0;
        while[smileIdx<count smileTable;
            smileRow:smileTable smileIdx;
            sabrIV:@[.sabr.impliedVol[forwardValue;smileRow`strike;smileRow`expiry;];sabrParams;{x}];
            if[10h=type sabrIV;
                failedCount+:1];
            if[not 10h=type sabrIV;
                modelVols:modelVols,sabrIV;
                marketVols:marketVols,smileRow`marketImpliedVolatility;
                weightVal:$[`weight in key smileRow;smileRow`weight;1f];
                weights:weights,weightVal];
            smileIdx+:1];
        pricedCount:count modelVols;
        rmseVal:$[pricedCount>0;.objective.rmse[modelVols;marketVols];0Nf];
        maeVal:$[pricedCount>0;.objective.mae[modelVols;marketVols];0Nf];
        diffs:modelVols-marketVols;
        sseVal:$[pricedCount>0;sum diffs*diffs;0Nf];
        wSseVal:$[pricedCount>0;sum weights*diffs*diffs;0Nf];
        calibrationTable:calibrationTable,enlist `calibrationId`alpha`beta`rho`nu`rmse`mae`totalSse`weightedSse`pricedRows`failedRows`status`errorMessage!(
            gridIdx+1;paramRow`alpha;paramRow`beta;paramRow`rho;paramRow`nu;
            rmseVal;maeVal;sseVal;wSseVal;pricedCount;failedCount;`OK;"");
        gridIdx+:1];
    calibrationTable
 };

.sabr.bestCalibration:{[calibrationTable]
    statusCol:calibrationTable`status;
    okMask:statusCol=`OK;
    okRows:calibrationTable where okMask;
    if[0=count okRows; '"No successful SABR calibration results"];
    rmseValues:okRows`rmse;
    bestIdx:rmseValues?min rmseValues;
    okRows bestIdx
 };

/ --- Pricing with SABR ---

.sabr.priceEuropeanWithSabr:{[optionRow;marketData;calibrationResult;configDict]
    spotVal:marketData`spot;
    riskFreeRate:marketData`riskFreeRate;
    dividendYield:marketData`dividendYield;
    expiryVal:optionRow`expiry;
    / Forward price
    rMinusQ:riskFreeRate-dividendYield;
    forwardValue:spotVal*exp rMinusQ*expiryVal;
    / SABR params from calibration
    sabrParams:`alpha`beta`rho`nu!(calibrationResult`alpha;calibrationResult`beta;calibrationResult`rho;calibrationResult`nu);
    / SABR implied vol
    sabrIV:@[.sabr.impliedVol[forwardValue;optionRow`strike;expiryVal;];sabrParams;{x}];
    if[10h=type sabrIV;
        :`tradeId`modelName`unitPrice`impliedVolatility`notionalPrice`status`statusMessage!(
            optionRow`optionId;`sabr;0Nf;0Nf;0Nf;`ERROR;sabrIV)];
    / BS price with SABR vol
    bsPrice:.validation.blackScholesClosedForm[optionRow`optionType;spotVal;optionRow`strike;expiryVal;riskFreeRate;dividendYield;sabrIV];
    notionalVal:$[`notional in key optionRow;optionRow`notional;1f];
    `tradeId`modelName`unitPrice`impliedVolatility`notionalPrice`status`statusMessage!(
        optionRow`optionId;`sabr;bsPrice;sabrIV;bsPrice*notionalVal;`OK;"")
 };
