\l core/init.q
/ shortVariance entry gate: open when IV > forecastVol + entryMargin, else FLAT.

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `SV_GATE;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;100;200;0f;300f;`linear;1b;1b);

openPath:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.30;7;1f%252f;0.02;0f;42);
closedPath:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(
    100f;0f;0.15;7;1f%252f;0.02;0f;42);

stratCfg:.strategy.defaultConfig `shortVariance;
stratCfg:@[stratCfg;(`forecastVol;`entryMargin;`stepYears);:;(0.20;0.02;1f%252f)];

openBundle:.strategy.runAndSummarize[`shortVariance;trade;openPath;bsModel;fdmCfg;stratCfg];
openSummary:openBundle`summary;
openResult:openBundle`result;
.testutil.assertTrue[openSummary[`gateOpen];"gate open when IV 0.30 > 0.20+0.02"];
.testutil.assertTrue[(openSummary`premiumCollected)>0f;"premium collected when gate open"];
.testutil.assertTrue[all (openResult`status)=`OK;"all rows OK when gate open"];

closedBundle:.strategy.runAndSummarize[`shortVariance;trade;closedPath;bsModel;fdmCfg;stratCfg];
closedSummary:closedBundle`summary;
closedResult:closedBundle`result;
.testutil.assertTrue[not closedSummary`gateOpen;"gate closed when IV 0.15 < 0.20+0.02"];
.testutil.assertTrue[0f=closedSummary`premiumCollected;"premium = 0 when gate closed"];
.testutil.assertTrue[0f=closedSummary`totalPnl;"totalPnl = 0 when gate closed"];
.testutil.assertTrue[`flat=closedSummary`status;"summary status flat when gate closed"];
.testutil.assertTrue[all (closedResult`status)=`flat;"every row flat when gate closed"];
.testutil.assertTrue[7=closedSummary`steps;"steps = path length even on flat path"];
.testutil.assertNear[closedSummary`impliedVolAtEntry;0.15;1e-12;"impliedVolAtEntry recorded even when gate closed"];

borderCfg:@[stratCfg;`forecastVol;:;0.25];
borderBundle:.strategy.runAndSummarize[`shortVariance;trade;closedPath;bsModel;fdmCfg;borderCfg];
.testutil.assertTrue[not borderBundle[`summary;`gateOpen];"IV 0.15 vs forecast+margin 0.27 -> closed"];

-1 "PASS test_strategy_short_variance_gate: openPremium=",string[openSummary`premiumCollected],", closedPnl=",string[closedSummary`totalPnl];
