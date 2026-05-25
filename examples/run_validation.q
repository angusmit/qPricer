/ =============================================================================
/ Example: Run validation against Black-Scholes closed form
/ =============================================================================
/ Validates FDM prices against analytical BS, checks put-call parity,
/ and runs a grid convergence test.
/ Run from the project root: q examples/run_validation.q
/ =============================================================================

\l lib/init.q

/ Define trades
callTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    10;`AAPL;`equityOption;`european;`call;100f;1f;1000000f
);

putTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    11;`AAPL;`equityOption;`european;`put;100f;1f;1000000f
);

/ Define flat market data
marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(
    `AAPL;100f;0.05;0f;0.2
);

/ Create model
model:.model.createBlackScholesModel[];

/ Configure solver
config:.config.createFiniteDifferenceConfig[`method`numberOfSpotSteps`numberOfTimeSteps`maximumSpot!(
    `explicit;200;1000;300f
)];

/ --- Validate call ---
-1 "=== European Call Validation ===";
callValidation:.validation.validateEuropeanOption[callTrade;marketData;model;config];
show callValidation;

/ --- Validate put ---
-1 "\n=== European Put Validation ===";
putValidation:.validation.validateEuropeanOption[putTrade;marketData;model;config];
show putValidation;

/ --- Put-call parity ---
-1 "\n=== Put-Call Parity Check ===";
parityResult:.validation.checkPutCallParity[callTrade;putTrade;marketData;model;config];
show parityResult;

/ --- Grid convergence ---
-1 "\n=== Grid Convergence Test ===";
configList:{
    .config.createFiniteDifferenceConfig[`method`numberOfSpotSteps`numberOfTimeSteps`maximumSpot!(
        `explicit;x;y;300f
    )]
 } .' flip (50 100 200 400;250 500 1000 2000);

convergenceResult:.validation.runGridConvergenceTest[callTrade;marketData;model;configList];
show convergenceResult;

-1 "\nValidation complete.";
