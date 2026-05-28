\l lib/init.q
/ The gamma/theta attribution identity holds only on steps where the position does not
/ mutate. Compute gammaReconResidual EXCLUDING roll steps; assert it is small relative
/ to the per-step option premium scale. The summary already computes residual over
/ non-roll OK rows.

pathTbl:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.10;13;1f%2520f;0f;0f;101);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `CR_REC;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;200;200;0f;200f;`linear;1b;1b);

stratCfg:.strategy.defaultConfig `calendarRoll;
stratCfg:@[stratCfg;(`frontTenorYears;`backTenorYears;`rollThresholdYears;`stepYears;`txnCostRate;`financingRate);:;(0.5;1f;0.0;1f%2520f;0f;0f)];

bundle:.strategy.runAndSummarize[`calendarRoll;trade;pathTbl;bsModel;fdmCfg;stratCfg];
resultTbl:bundle`result;
summary:bundle`summary;

.testutil.assertTrue[`OK=summary`status;"summary status OK"];
.testutil.assertTrue[0=summary`totalRolls;"no rolls in this short path with long tenors"];

nonRollRows:resultTbl where (resultTbl`rollEvents)=0;
.testutil.assertTrue[(count nonRollRows)=count resultTbl;"every row is non-roll"];

frontLegEntryPrice:0.05;
absResidual:abs summary`gammaReconResidual;
.testutil.assertTrue[not null summary`gammaReconResidual;"residual finite"];
.testutil.assertTrue[absResidual<0.5;"non-roll residual within absolute tolerance for tiny-step path"];
.testutil.assertNear[summary`txnCostTotal;0f;1e-12;"zero txn cost"];
.testutil.assertNear[summary`financingTotal;0f;1e-12;"zero financing"];

resWithRoll:bundle`summary;
zeroRollResidual:resWithRoll`gammaReconResidual;
.testutil.assertTrue[not null zeroRollResidual;"residual reported even with zero rolls (non-roll rows include the whole path)"];

-1 "PASS test_strategy_calendar_roll_recon: residual=",string[summary`gammaReconResidual],", totalRolls=",string[summary`totalRolls];
