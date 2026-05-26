/ =============================================================================
/ Test: European put pricing
/ =============================================================================
/ Verifies that the FDM put price is within tolerance of BS closed form.
/ =============================================================================

\l lib/init.q

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    101;`TEST;`equityOption;`european;`put;100f;1f;1f
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
bsPrice:.validation.blackScholesClosedForm[`put;100f;100f;1f;0.05;0f;0.2];

/ Check tolerance
tolerance:0.25;
absError:abs fdmPrice - bsPrice;

$[absError < tolerance;
    -1 "PASS test_european_put: FDM=",string[fdmPrice],", BS=",string[bsPrice],", error=",string absError;
    [
        -2 "FAIL test_european_put: FDM=",string[fdmPrice],", BS=",string[bsPrice],", error=",string[absError]," > tolerance=",string tolerance;
        '"test_european_put FAILED"
    ]
];
