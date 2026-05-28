\l lib/init.q
/ Drive gammaScalp through the Barchart path adapter using a small synthetic
/ normalised table. Structural-only assertion (the time-to-expiry doesn't decay
/ across snapshots in v0.42; this is an integration test, not a realistic
/ historical replay).

synthNormalised:([] snapshotDate:2024.01.02 2024.01.03 2024.01.04 2024.01.05 2024.01.08;
                    underlying:5#`aapl;
                    contractId:5#`aapl_C100;
                    expiryDate:5#2024.02.16;
                    optionType:5#`call;
                    strike:5#100f;
                    spot:185.0 186.5 184.2 187.1 186.0;
                    marketPrice:5.5 6.1 5.4 6.2 5.9;
                    impliedVolatility:0.28 0.27 0.29 0.275 0.28;
                    status:5#`OK);

pathTbl:.strategy.path.fromBarchart[synthNormalised;`aapl_C100];

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    `T_BC;`aapl;`equityOption;`european;`call;100f;0.125;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `crankNicolson;100;150;0f;400f;`linear;1b;1b);

stratCfg:.strategy.defaultConfig `gammaScalp;
stratCfg:@[stratCfg;(`rebalanceMode;`rebalanceInterval;`stepYears);:;(`interval;1;1f%252f)];

bundle:.strategy.runAndSummarize[`gammaScalp;trade;pathTbl;bsModel;fdmCfg;stratCfg];
resultTbl:bundle`result;
summary:bundle`summary;

.testutil.assertTrue[5=count resultTbl;"result rows match path rows"];
.testutil.assertTrue[all (resultTbl`status)=`OK;"every OK"];
.testutil.assertTrue[(summary[`status]) in `OK`warning`ERROR;"summary status is a recognised symbol"];
.testutil.assertTrue[not null summary`totalPnl;"totalPnl finite"];
.testutil.assertTrue[(summary`steps)=5;"summary steps = 5"];

-1 "PASS test_strategy_gamma_scalp_barchart: steps=",string[summary`steps],", totalPnl=",string[summary`totalPnl];
