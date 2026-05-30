/ validation.q — Black-Scholes closed form, Greeks, and validation tools
/ Normal CDF: Abramowitz & Stegun 26.2.17, accuracy < 7.5e-8

/ --- Normal distribution ---

.validation.__normalCdf:{[x]
    if[x<0f; :1f-.validation.__normalCdf neg x];
    t:1f%(1f+0.2316419*x);
    pdf:exp[neg 0.5*x*x]%sqrt 2f*acos neg 1f;
    t2:t*t; t3:t2*t; t4:t3*t; t5:t4*t;
    1f-pdf*(0.31938153*t)+(-0.356563782*t2)+(1.781477937*t3)+(-1.821255978*t4)+1.330274429*t5
 };

.validation.__normalPdf:{[x]
    exp[neg 0.5*x*x]%sqrt 2f*acos neg 1f
 };

/ --- Black-Scholes helpers ---

.validation.__d1:{[spot;strike;expiry;rate;divY;vol]
    sqrtT:sqrt expiry;
    (log[spot%strike]+((rate-divY)+0.5*vol*vol)*expiry)%(vol*sqrtT)
 };

.validation.__d2:{[d1Val;vol;expiry] d1Val-vol*sqrt expiry};

/ --- Closed-form price ---

.validation.blackScholesClosedForm:{[optType;spot;strike;expiry;rate;divY;vol]
    d1Val:.validation.__d1[spot;strike;expiry;rate;divY;vol];
    d2Val:.validation.__d2[d1Val;vol;expiry];
    rDisc:exp neg rate*expiry;
    qDisc:exp neg divY*expiry;
    isCall:optType~`call;
    if[isCall; :(spot*qDisc*.validation.__normalCdf d1Val)-strike*rDisc*.validation.__normalCdf d2Val];
    (strike*rDisc*.validation.__normalCdf neg d2Val)-spot*qDisc*.validation.__normalCdf neg d1Val
 };

/ --- Analytical Greeks ---
/ Conventions: delta per 1 spot, gamma per 1 spot^2,
/ theta annual, vega per 1.00 vol, rho per 1.00 rate

.validation.blackScholesGreeks:{[optType;spot;strike;expiry;rate;divY;vol]
    d1Val:.validation.__d1[spot;strike;expiry;rate;divY;vol];
    d2Val:.validation.__d2[d1Val;vol;expiry];
    sqrtT:sqrt expiry;
    rDisc:exp neg rate*expiry;
    qDisc:exp neg divY*expiry;
    nd1:.validation.__normalCdf d1Val;
    nd2:.validation.__normalCdf d2Val;
    pdf1:.validation.__normalPdf d1Val;
    isCall:optType~`call;
    / Delta
    delta:$[isCall; qDisc*nd1; neg qDisc*.validation.__normalCdf neg d1Val];
    / Gamma (same for call and put)
    gamma:qDisc*pdf1%(spot*vol*sqrtT);
    / Theta (annual)
    thetaTerm1:neg spot*qDisc*pdf1*vol%(2f*sqrtT);
    callTheta:thetaTerm1 - (rate*strike*rDisc*nd2) - divY*spot*qDisc*nd1;
    putTheta:thetaTerm1 + (rate*strike*rDisc*.validation.__normalCdf neg d2Val) - divY*spot*qDisc*.validation.__normalCdf neg d1Val;
    theta:$[isCall;callTheta;putTheta];
    / Vega (per 1.00 absolute vol)
    vega:spot*qDisc*pdf1*sqrtT;
    / Rho (per 1.00 absolute rate)
    rho:$[isCall; strike*expiry*rDisc*nd2; neg strike*expiry*rDisc*.validation.__normalCdf neg d2Val];
    `delta`gamma`theta`vega`rho!(delta;gamma;theta;vega;rho)
 };

/ --- Validation: price ---

.validation.validateEuropeanOption:{[trade;marketData;model;config]
    fdm:(.engine.priceOption[trade;marketData;model;config])`unitPrice;
    bs:.validation.blackScholesClosedForm[trade`optionType;marketData`spot;trade`strike;
        trade`expiry;marketData`riskFreeRate;marketData`dividendYield;marketData`volatility];
    ([] tradeId:enlist trade`tradeId; optionType:enlist trade`optionType;
        fdmUnitPrice:enlist fdm; closedFormUnitPrice:enlist bs;
        absoluteError:enlist abs fdm-bs; relativeError:enlist .utilities.relativeError[fdm;bs])
 };

/ --- Validation: Greeks ---

.validation.validateGreeks:{[trade;marketData;model;config]
    fdmGreeks:.greeks.calculateGreeks[trade;marketData;model;config];
    bsGreeks:.validation.blackScholesGreeks[trade`optionType;marketData`spot;trade`strike;
        trade`expiry;marketData`riskFreeRate;marketData`dividendYield;marketData`volatility];
    greekNames:`delta`gamma`theta`vega`rho;
    fdmVals:(first each fdmGreeks greekNames);
    bsVals:bsGreeks greekNames;
    absErrs:abs fdmVals - bsVals;
    relErrs:{.utilities.relativeError[x;y]}'[fdmVals;bsVals];
    ([] tradeId:(count greekNames)#trade`tradeId;
        optionType:(count greekNames)#trade`optionType;
        greek:greekNames;
        fdmValue:fdmVals;
        closedFormValue:bsVals;
        absoluteError:absErrs;
        relativeError:relErrs)
 };

/ --- Put-call parity ---

.validation.checkPutCallParity:{[callTrade;putTrade;marketData;model;config]
    cP:(.engine.priceOption[callTrade;marketData;model;config])`unitPrice;
    pP:(.engine.priceOption[putTrade;marketData;model;config])`unitPrice;
    actual:cP-pP;
    rDisc:exp neg marketData[`riskFreeRate]*callTrade`expiry;
    qDisc:exp neg marketData[`dividendYield]*callTrade`expiry;
    theoretical:(marketData[`spot]*qDisc)-callTrade[`strike]*rDisc;
    `callPrice`putPrice`actualDifference`theoreticalDifference`absoluteError!(
        cP;pP;actual;theoretical;abs actual-theoretical)
 };

/ --- Grid convergence ---

.validation.runGridConvergenceTest:{[trade;marketData;model;configList]
    bs:.validation.blackScholesClosedForm[trade`optionType;marketData`spot;trade`strike;
        trade`expiry;marketData`riskFreeRate;marketData`dividendYield;marketData`volatility];
    results:{[trade;marketData;model;bs;cfg]
        fdm:(.engine.priceOption[trade;marketData;model;cfg])`unitPrice;
        `numberOfSpotSteps`numberOfTimeSteps`fdmPrice`closedFormPrice`absoluteError`relativeError!(
            cfg`numberOfSpotSteps; cfg`numberOfTimeSteps; fdm; bs;
            abs fdm-bs; .utilities.relativeError[fdm;bs])
    }[trade;marketData;model;bs] each configList;
    columns:`numberOfSpotSteps`numberOfTimeSteps`fdmPrice`closedFormPrice`absoluteError`relativeError;
    flip columns!flip value each results
 };
