/ test_portfolio_heston_products.q
\l core/init.q
hestonParams:`initialVariance`longRunVariance`meanReversion`volOfVol`correlation`riskFreeRate`dividendYield!(0.04;0.04;2.0;0.3;-0.7;0.05;0.0);
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`momentMatching`confidenceLevel!(10000;50;42;0b;0b;0.95);
configDict:`mcConfig`hestonParams`modelType!(mcConfig;hestonParams;`heston);

/ European call under Heston
callTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;100000f);
callResult:.heston.priceEuropean[callTrade;mkt;configDict];
.testutil.assertTrue[callResult[`status]~`OK;"Heston call OK"];
.testutil.assertTrue[callResult[`unitPrice]>0f;"Heston call price positive"];
.testutil.assertNear[callResult`notionalPrice;callResult[`unitPrice]*100000f;1f;"notional correct"];

/ American under Heston should error via portfolio routing
americanTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate`averageType`averagingStyle`observationCount`basketSymbols`basketWeights`lookbackStyle!(
    2;`AAPL;`equityOption;`american;`call;100f;1f;100000f;`none;0Nf;0f;`none;`none;0N;`symbol$();`float$();`none);
mdl:.model.createBlackScholesModel[];
cfg:.config.defaultPricingConfig[];
hestonCfg:cfg,`mcConfig`hestonParams`modelType!(mcConfig;hestonParams;`heston);
tradeTable:enlist americanTrade;
portfolioResult:.portfolio.priceTrades[tradeTable;mkt;mdl;hestonCfg];
americanRow:portfolioResult 0;
.testutil.assertTrue[americanRow[`status]~`ERROR;"American Heston errors"];

-1 "PASS test_portfolio_heston_products: hestonCallPrice=",string callResult`unitPrice;
