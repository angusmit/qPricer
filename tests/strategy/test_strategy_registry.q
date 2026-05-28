\l lib/init.q
/ The generic engine .strategy.run must dispatch through the registry without any
/ strategy-specific branching. Demonstrate by registering a trivial dummy strategy
/ and running it through the SAME driver that runs gammaScalp.

dummyInit:{[trade;firstStep;model;fdmConfig;stratCfg]
    initialVal:`float$stratCfg`startValue;
    rowEmit:`stepIndex`stepDate`spot`value`status`message!(
        firstStep`stepIndex;firstStep`stepDate;firstStep`spot;initialVal;`OK;"");
    `cumulativeValue`rowEmit!(initialVal;rowEmit)
 };

dummyStep:{[state;marketStep;trade;model;fdmConfig;stratCfg]
    incrementVal:`float$stratCfg`increment;
    newValue:state[`cumulativeValue]+incrementVal;
    rowEmit:`stepIndex`stepDate`spot`value`status`message!(
        marketStep`stepIndex;marketStep`stepDate;marketStep`spot;newValue;`OK;"");
    @[state;(`cumulativeValue;`rowEmit);:;(newValue;rowEmit)]
 };

dummySummary:{[resultTable;stratCfg]
    `strategyName`steps`totalValue!(`dummy;count resultTable;sum resultTable`value)
 };

dummyDefaultCfg:{[] `startValue`increment!(0f;1f)};

.strategy.register[`dummy;dummyInit;dummyStep;dummySummary;dummyDefaultCfg];

.testutil.assertTrue[`dummy in .strategy.registeredStrategies[];"dummy registered"];
.testutil.assertTrue[`gammaScalp in .strategy.registeredStrategies[];"gammaScalp still registered"];
.testutil.assertTrue[(.strategy.defaultConfig`dummy)~`startValue`increment!(0f;1f);"defaultConfig returns dummy keys"];

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(`T1;`X;`equityOption;`european;`call;100f;0.25;1f);
pathTbl:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.20;5;1f%252f;0f;0f;7);

emptyDict:()!();
resultTbl:.strategy.run[`dummy;trade;pathTbl;emptyDict;emptyDict;.strategy.defaultConfig`dummy];
.testutil.assertTrue[5=count resultTbl;"dummy run produced 5 rows"];
.testutil.assertTrue[(`float$til 5)~resultTbl`value;"dummy step counter matches expected sequence"];

bundle:.strategy.runAndSummarize[`dummy;trade;pathTbl;emptyDict;emptyDict;.strategy.defaultConfig`dummy];
.testutil.assertTrue[`result`summary~asc key bundle;"runAndSummarize returns result+summary"];
.testutil.assertTrue[`dummy=bundle[`summary;`strategyName];"summary reports dummy strategyName"];

unknownResult:.[.strategy.run;(`notARegisteredStrategy;trade;pathTbl;emptyDict;emptyDict;emptyDict);{`ERROR}];
.testutil.assertTrue[unknownResult~`ERROR;"unknown strategyName rejected"];

-1 "PASS test_strategy_registry";
