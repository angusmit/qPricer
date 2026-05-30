\l core/init.q
/ Construct an explicit path that ramps spot above the up-and-out barrier
/ to guarantee a knock event mid-run.
nSteps:5;
spots:100 104 108 112 116f;
pathTbl:flip `stepIndex`stepDate`spot`volatility`riskFreeRate`dividendYield`marketPrice`status!(
    til nSteps;
    2024.01.01+til nSteps;
    spots;
    nSteps#0.20;
    nSteps#0.02;
    nSteps#0f;
    nSteps#0Nf;
    nSteps#`OK);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel!(
    `BH_K;`X;`equityOption;`european;`call;100f;0.25;1f;`upAndOut;110f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;30;500;0f;300f;`linear;1b;0b);
stratCfg:@[.strategy.defaultConfig `barrierHedge;`stepYears;:;1f%252f];
bundle:.strategy.runAndSummarize[`barrierHedge;trade;pathTbl;bsModel;fdmCfg;stratCfg];
resultTbl:bundle`result;
summary:bundle`summary;
.testutil.assertTrue[summary`knockedOut;"barrier breached on aggressive path"];
knockedRows:resultTbl where resultTbl`knockedOut;
.testutil.assertTrue[0<count knockedRows;"at least one knocked row"];
postKnockRows:resultTbl where (resultTbl`stepIndex)>summary`knockStep;
.testutil.assertTrue[all 0f=postKnockRows`optionPrice;"option price 0 past knock (rebate)"];
.testutil.assertTrue[all 0f=postKnockRows`positionValue;"positionValue 0 past knock"];
.testutil.assertTrue[all 0f=postKnockRows`hedgePosition;"hedge closed past knock"];
.testutil.assertTrue[all 0f=postKnockRows`hedgeTrade;"no hedge trades past knock"];
.testutil.assertTrue[all (postKnockRows`status)=`OK;"post-knock rows OK"];
-1 "PASS test_strategy_barrier_hedge_knockout: knockStep=",string[summary`knockStep],", barrierLevel=",string[summary`barrierLevel];
