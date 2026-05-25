/ =============================================================================
/ Example: Calculate Greeks
/ =============================================================================
/ Computes delta, gamma, theta, vega, and rho for a European call.
/ Run from the project root: q examples/calculate_greeks.q
/ =============================================================================

\l lib/init.q

/ Define the trade
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    4;
    `AAPL;
    `equityOption;
    `european;
    `call;
    100f;
    1f;
    1000000f
 );

/ Define flat market data
marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(
    `AAPL;
    100f;
    0.05;
    0f;
    0.2
 );

/ Create model
model:.model.createBlackScholesModel[];

/ Configure solver (moderate grid for reasonable Greek accuracy)
config:.config.createFiniteDifferenceConfig[`method`numberOfSpotSteps`numberOfTimeSteps`maximumSpot!(
    `explicit;200;1000;300f
 )];

/ Calculate Greeks
-1 "Calculating Greeks for European call...";
greeks:.greeks.calculateGreeks[trade;marketData;model;config];
show greeks;

-1 "\nGreek values:";
-1 "  Delta: ",string first greeks`delta;
-1 "  Gamma: ",string first greeks`gamma;
-1 "  Theta: ",string first greeks`theta;
-1 "  Vega:  ",string first greeks`vega;
-1 "  Rho:   ",string first greeks`rho;
