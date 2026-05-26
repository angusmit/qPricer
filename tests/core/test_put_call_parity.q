/ =============================================================================
/ Test: Put-call parity
/ =============================================================================
/ Verifies that C - P ≈ S*exp(-q*T) - K*exp(-r*T)
/ =============================================================================

\l lib/init.q

callTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    102;`TEST;`equityOption;`european;`call;100f;1f;1f
 );

putTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    103;`TEST;`equityOption;`european;`put;100f;1f;1f
 );

marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(
    `TEST;100f;0.05;0f;0.2
 );

model:.model.createBlackScholesModel[];

config:.config.createFiniteDifferenceConfig[`method`numberOfSpotSteps`numberOfTimeSteps`maximumSpot!(
    `explicit;200;2000;300f
 )];

/ Check parity
parityResult:.validation.checkPutCallParity[callTrade;putTrade;marketData;model;config];

tolerance:0.5;
absError:parityResult`absoluteError;

$[absError < tolerance;
    -1 "PASS test_put_call_parity: actual=",string[parityResult`actualDifference],", theoretical=",string[parityResult`theoreticalDifference],", error=",string absError;
    [
        -2 "FAIL test_put_call_parity: error=",string[absError]," > tolerance=",string tolerance;
        '"test_put_call_parity FAILED"
    ]
];
