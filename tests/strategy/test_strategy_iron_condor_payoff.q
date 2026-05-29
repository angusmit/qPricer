\l lib/init.q
/ Analytic identity check: maxProfit/maxLoss/breakevens against independent
/ formula based on spread widths and netCredit.
pathTbl:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.25;5;1f%252f;0.02;0f;7);
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `IC_P;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;80;150;0f;300f;`linear;1b;1b);
stratCfg:@[.strategy.defaultConfig `ironCondor;`stepYears;:;1f%252f];
bundle:.strategy.runAndSummarize[`ironCondor;trade;pathTbl;bsModel;fdmCfg;stratCfg];
summary:bundle`summary;
spot0:100f;
shortPutStrike:spot0*1f-stratCfg`shortPutPct;
longPutStrike:spot0*1f-stratCfg`longPutPct;
shortCallStrike:spot0*1f+stratCfg`shortCallPct;
longCallStrike:spot0*1f+stratCfg`longCallPct;
expectedPutWidth:shortPutStrike-longPutStrike;
expectedCallWidth:longCallStrike-shortCallStrike;
expectedWidestWidth:expectedPutWidth|expectedCallWidth;
expectedMaxProfit:summary`netCredit;
expectedMaxLoss:expectedWidestWidth-summary`netCredit;
expectedLowerBE:shortPutStrike-summary`netCredit;
expectedUpperBE:shortCallStrike+summary`netCredit;
.testutil.assertTrue[1e-12>abs (summary`putSpreadWidth)-expectedPutWidth;"putSpreadWidth analytic"];
.testutil.assertTrue[1e-12>abs (summary`callSpreadWidth)-expectedCallWidth;"callSpreadWidth analytic"];
.testutil.assertTrue[1e-12>abs (summary`widestSpreadWidth)-expectedWidestWidth;"widestSpreadWidth analytic"];
.testutil.assertTrue[1e-12>abs (summary`maxProfit)-expectedMaxProfit;"maxProfit analytic"];
.testutil.assertTrue[1e-12>abs (summary`maxLoss)-expectedMaxLoss;"maxLoss analytic"];
.testutil.assertTrue[1e-12>abs (summary`lowerBreakeven)-expectedLowerBE;"lowerBreakeven analytic"];
.testutil.assertTrue[1e-12>abs (summary`upperBreakeven)-expectedUpperBE;"upperBreakeven analytic"];
-1 "PASS test_strategy_iron_condor_payoff: lowerBE=",string[summary`lowerBreakeven],", upperBE=",string[summary`upperBreakeven];
