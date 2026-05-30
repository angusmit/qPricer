\l core/init.q
/ Variance risk premium: when realized vol falls well below implied, short straddle
/ should profit (positive totalPnl). When realized vol exceeds implied, totalPnl
/ should fall. We synthesise two paths with the same implied vol but different
/ realized vols and compare summaries.

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `SV_PREM;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;100;200;0f;300f;`linear;1b;1b);
stratCfg:.strategy.defaultConfig `shortVariance;
stratCfg:@[stratCfg;(`forecastVol;`entryMargin;`stepYears);:;(0.10;0.02;1f%252f)];

quietPath:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed`startDate!(
    100f;0f;0.05;11;1f%252f;0.02;0f;7;2024.01.01);
quietPath:update volatility:count[quietPath]#0.30 from quietPath;
noisyPath:.strategy.path.fromSynthetic `spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed`startDate!(
    100f;0f;0.50;11;1f%252f;0.02;0f;7;2024.01.01);
noisyPath:update volatility:count[noisyPath]#0.30 from noisyPath;

quietBundle:.strategy.runAndSummarize[`shortVariance;trade;quietPath;bsModel;fdmCfg;stratCfg];
noisyBundle:.strategy.runAndSummarize[`shortVariance;trade;noisyPath;bsModel;fdmCfg;stratCfg];
quietSummary:quietBundle`summary;
noisySummary:noisyBundle`summary;

.testutil.assertTrue[quietSummary`gateOpen;"quiet-path gate open (IV 0.30 > 0.10+0.02)"];
.testutil.assertTrue[noisySummary`gateOpen;"noisy-path gate open (same IV)"];

.testutil.assertNear[quietSummary`impliedVolAtEntry;0.30;1e-12;"quiet implied 0.30"];
.testutil.assertNear[noisySummary`impliedVolAtEntry;0.30;1e-12;"noisy implied 0.30"];
.testutil.assertTrue[(quietSummary`realizedVol)<noisySummary`realizedVol;"realized vol lower on quiet path"];
.testutil.assertTrue[(quietSummary`varianceRiskPremium)>noisySummary`varianceRiskPremium;"VRP (impl - real) larger on quiet path"];

.testutil.assertTrue[(quietSummary`totalPnl)>noisySummary`totalPnl;"quiet path totalPnl > noisy path (short variance harvested)"];

-1 "PASS test_strategy_short_variance_premium: quietPnl=",string[quietSummary`totalPnl],", noisyPnl=",string[noisySummary`totalPnl],", quietVRP=",string[quietSummary`varianceRiskPremium],", noisyVRP=",string[noisySummary`varianceRiskPremium];
