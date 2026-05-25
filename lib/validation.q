/ =============================================================================
/ validation.q — Validation against Black-Scholes closed form
/ =============================================================================
/ Provides Black-Scholes closed-form pricing for European options and
/ validation tools to compare FDM results against analytical solutions.
/ All functions live in the .validation namespace.
/ =============================================================================

/ -----------------------------------------------------------------------------
/ Public functions — closed-form pricing
/ -----------------------------------------------------------------------------

/ Black-Scholes closed-form price for a European option with dividends.
/ Call = S*exp(-q*T)*N(d1) - K*exp(-r*T)*N(d2)
/ Put  = K*exp(-r*T)*N(-d2) - S*exp(-q*T)*N(-d1)
/ Parameters:
/   optionType    - `call or `put
/   spot          - current spot price
/   strike        - strike price
/   expiry        - time to expiry in years
/   riskFreeRate  - annualised risk-free rate
/   dividendYield - annualised continuous dividend yield
/   volatility    - annualised volatility
/ Returns: option price (float)
.validation.blackScholesClosedForm:{[optionType;spot;strike;expiry;riskFreeRate;dividendYield;volatility]
    d1:.validation.__blackScholesD1[spot;strike;expiry;riskFreeRate;dividendYield;volatility];
    d2:.validation.__blackScholesD2[d1;volatility;expiry];

    rateDiscount:exp neg riskFreeRate * expiry;
    dividendDiscount:exp neg dividendYield * expiry;

    $[optionType ~ `call;
        (spot * dividendDiscount * .validation.__normalCdf[d1]) - strike * rateDiscount * .validation.__normalCdf[d2];
      optionType ~ `put;
        (strike * rateDiscount * .validation.__normalCdf[neg d2]) - spot * dividendDiscount * .validation.__normalCdf[neg d1];
      '"Unsupported optionType for closed form: ",string optionType
    ]
 };

/ -----------------------------------------------------------------------------
/ Public functions — validation
/ -----------------------------------------------------------------------------

/ Validate FDM price against Black-Scholes closed form for a European option.
/ Parameters:
/   trade      - trade dictionary
/   marketData - market data dictionary
/   model      - model dictionary
/   config     - finite-difference config dictionary
/ Returns: one-row table with fdmUnitPrice, closedFormUnitPrice,
/          absoluteError, relativeError

.validation.validateEuropeanOption:{[trade;marketData;model;config]
    / FDM price
    priceResult:.engine.priceOption[trade;marketData;model;config];
    fdmUnitPrice:priceResult`unitPrice;

    / Closed-form price
    closedFormUnitPrice:.validation.blackScholesClosedForm[
        trade`optionType;
        marketData`spot;
        trade`strike;
        trade`expiry;
        marketData`riskFreeRate;
        marketData`dividendYield;
        marketData`volatility
    ];

    absError:abs fdmUnitPrice - closedFormUnitPrice;
    relError:.utilities.relativeError[fdmUnitPrice;closedFormUnitPrice];

    ([]
        tradeId:enlist trade`tradeId;
        optionType:enlist trade`optionType;
        fdmUnitPrice:enlist fdmUnitPrice;
        closedFormUnitPrice:enlist closedFormUnitPrice;
        absoluteError:enlist absError;
        relativeError:enlist relError
    )
 };

/ Check put-call parity: C - P = S*exp(-q*T) - K*exp(-r*T)
/ Parameters:
/   callTrade  - trade dictionary for the call
/   putTrade   - trade dictionary for the put
/   marketData - market data dictionary
/   model      - model dictionary
/   config     - finite-difference config dictionary
/ Returns: dictionary with actual difference, theoretical difference,
/          absolute error
.validation.checkPutCallParity:{[callTrade;putTrade;marketData;model;config]
    callResult:.engine.priceOption[callTrade;marketData;model;config];
    putResult:.engine.priceOption[putTrade;marketData;model;config];

    callPrice:callResult`unitPrice;
    putPrice:putResult`unitPrice;
    actualDifference:callPrice - putPrice;

    spot:marketData`spot;
    strike:callTrade`strike;
    expiry:callTrade`expiry;
    riskFreeRate:marketData`riskFreeRate;
    dividendYield:marketData`dividendYield;

    rateDiscount:exp neg riskFreeRate * expiry;
    dividendDiscount:exp neg dividendYield * expiry;
    theoreticalDifference:(spot * dividendDiscount) - strike * rateDiscount;

    absError:abs actualDifference - theoreticalDifference;

    `callPrice`putPrice`actualDifference`theoreticalDifference`absoluteError!(
        callPrice;
        putPrice;
        actualDifference;
        theoreticalDifference;
        absError
    )
 };

/ Run a grid convergence test across multiple configurations.
/ Shows how FDM price converges to the closed-form as grid resolution increases.
/ Parameters:
/   trade      - trade dictionary
/   marketData - market data dictionary
/   model      - model dictionary
/   configList - list of config dictionaries with increasing resolution
/ Returns: table with numberOfSpotSteps, numberOfTimeSteps, fdmPrice,
/          closedFormPrice, absoluteError, relativeError
.validation.runGridConvergenceTest:{[trade;marketData;model;configList]
    / Closed-form price (constant across all configs)
    closedFormPrice:.validation.blackScholesClosedForm[
        trade`optionType;
        marketData`spot;
        trade`strike;
        trade`expiry;
        marketData`riskFreeRate;
        marketData`dividendYield;
        marketData`volatility
    ];

    / Run FDM for each config
    results:{[trade;marketData;model;closedFormPrice;config]
        priceResult:.engine.priceOption[trade;marketData;model;config];
        fdmPrice:priceResult`unitPrice;
        absError:abs fdmPrice - closedFormPrice;
        relError:.utilities.relativeError[fdmPrice;closedFormPrice];
        `numberOfSpotSteps`numberOfTimeSteps`fdmPrice`closedFormPrice`absoluteError`relativeError!(
            config`numberOfSpotSteps;
            config`numberOfTimeSteps;
            fdmPrice;
            closedFormPrice;
            absError;
            relError
        )
    }[trade;marketData;model;closedFormPrice] each configList;

    / Convert list of dictionaries to table
    resultKeys:`numberOfSpotSteps`numberOfTimeSteps`fdmPrice`closedFormPrice`absoluteError`relativeError;
    flip resultKeys!flip value each results
 };

/ -----------------------------------------------------------------------------
/ Internal functions — normal CDF and Black-Scholes helpers
/ -----------------------------------------------------------------------------

/ Approximate the standard normal cumulative distribution function N(x).
/ Uses the Abramowitz and Stegun approximation (equation 26.2.17).
/ Accuracy: absolute error < 7.5e-8.
/ Parameters:
/   x - input value
/ Returns: N(x) (float between 0 and 1)
.validation.__normalCdf:{[x]
    / Handle negative x by symmetry: N(x) = 1 - N(-x)
    if[x < 0f; :1f - .validation.__normalCdf[neg x]];

    / Polynomial coefficients (Abramowitz and Stegun 26.2.17)
    a1:0.319381530;
    a2:-0.356563782;
    a3:1.781477937;
    a4:-1.821255978;
    a5:1.330274429;
    p:0.2316419;

    / Transform variable
    t:1f % (1f + p * x);

    / Standard normal PDF: n(x) = exp(-x^2/2) / sqrt(2*pi)
    twoPi:2f * acos neg 1f;
    pdf:exp[neg 0.5 * x * x] % sqrt twoPi;

    / Polynomial in t (powers computed explicitly for clarity)
    t2:t * t;
    t3:t2 * t;
    t4:t3 * t;
    t5:t4 * t;
    poly:(a1 * t) + (a2 * t2) + (a3 * t3) + (a4 * t4) + a5 * t5;

    / N(x) = 1 - n(x) * poly(t)
    1f - pdf * poly
 };

/ Calculate Black-Scholes d1 parameter.
/ d1 = (log(S/K) + (r - q + 0.5*vol^2)*T) / (vol * sqrt(T))
/ Parameters:
/   spot, strike, expiry, riskFreeRate, dividendYield, volatility
/ Returns: d1 (float)
.validation.__blackScholesD1:{[spot;strike;expiry;riskFreeRate;dividendYield;volatility]
    sqrtExpiry:sqrt expiry;
    volSqrt:volatility * sqrtExpiry;
    logMoneyness:log spot % strike;
    driftTerm:((riskFreeRate - dividendYield) + 0.5 * volatility * volatility) * expiry;
    (logMoneyness + driftTerm) % volSqrt
 };

/ Calculate Black-Scholes d2 parameter.
/ d2 = d1 - vol * sqrt(T)
/ Parameters:
/   d1         - d1 value from __blackScholesD1
/   volatility - annualised volatility
/   expiry     - time to expiry in years
/ Returns: d2 (float)
.validation.__blackScholesD2:{[d1;volatility;expiry]
    d1 - volatility * sqrt expiry
 };

