/ =============================================================================
/ Test: Grid convergence
/ =============================================================================
/ Verifies that the FDM error decreases as grid resolution increases.
/ Uses a coarse and fine grid; fine grid error should be smaller.
/ =============================================================================

\l lib/init.q

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    104;`TEST;`equityOption;`european;`call;100f;1f;1f
 );

marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(
    `TEST;100f;0.05;0f;0.2
 );

model:.model.createBlackScholesModel[];

/ Coarse config
configCoarse:.config.createFiniteDifferenceConfig[`method`numberOfSpotSteps`numberOfTimeSteps`maximumSpot!(
    `explicit;50;250;300f
 )];

/ Fine config
configFine:.config.createFiniteDifferenceConfig[`method`numberOfSpotSteps`numberOfTimeSteps`maximumSpot!(
    `explicit;200;2000;300f
 )];

/ Run convergence test
configList:(configCoarse;configFine);
convergenceResult:.validation.runGridConvergenceTest[trade;marketData;model;configList];

-1 "Grid convergence results:";
show convergenceResult;

/ Check that fine grid has smaller error than coarse grid
errors:convergenceResult`absoluteError;
coarseError:errors 0;
fineError:errors 1;

$[fineError < coarseError;
    -1 "PASS test_grid_convergence: coarse error=",string[coarseError],", fine error=",string fineError;
    [
        -2 "FAIL test_grid_convergence: fine error (",string[fineError],") >= coarse error (",string[coarseError],")";
        '"test_grid_convergence FAILED"
    ]
];
