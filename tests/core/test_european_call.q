/ =============================================================================
/ Test: European call pricing
/ =============================================================================
/ Verifies that the FDM call price is within tolerance of BS closed form.
/ =============================================================================

\l lib/init.q

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    100;`TEST;`equityOption;`european;`call;100f;1f;1f
 );

marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(
    `TEST;100f;0.05;0f;0.2
 );

model:.model.createBlackScholesModel[];

config:.config.createFiniteDifferenceConfig[`method`numberOfSpotSteps`numberOfTimeSteps`maximumSpot!(
    `explicit;200;2000;300f
 )];

/ FDM price
priceResult:.engine.priceOption[trade;marketData;model;config];
fdmPrice:priceResult`unitPrice;

/ Closed-form price
bsPrice:.validation.blackScholesClosedForm[`call;100f;100f;1f;0.05;0f;0.2];

/ Check tolerance
tolerance:0.25;
absError:abs fdmPrice - bsPrice;

$[absError < tolerance;
    -1 "PASS test_european_call: FDM=",string[fdmPrice],", BS=",string[bsPrice],", error=",string absError;
    [
        -2 "FAIL test_european_call: FDM=",string[fdmPrice],", BS=",string[bsPrice],", error=",string[absError]," > tolerance=",string tolerance;
        '"test_european_call FAILED"
    ]
];
