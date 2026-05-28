/ test_timing.q - timing utilities
\l lib/init.q

/ 1. timeFunction with a simple function
timedResult:.timing.timeFunction[{[argVal] argVal*argVal};enlist 5];
.testutil.assertEqual[timedResult`result;25;"timeFunction result"];
.testutil.assertTrue[timedResult[`elapsedMs]>=0;"timeFunction elapsed non-negative"];

/ 2. timeFunction with engine pricing
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1f);
marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
bsModel:.model.createBlackScholesModel[];
fdmConfig:.config.defaultPricingConfig[];

timedPrice:.timing.timeFunction[.engine.priceOption;(trade;marketData;bsModel;fdmConfig)];
.testutil.assertTrue[99h=type timedPrice`result;"timed pricing returns dict"];
.testutil.assertTrue[timedPrice[`elapsedMs]>0;"pricing takes measurable time"];

/ 3. benchmarkPricingGrid
configList:(
    .config.mergeConfig[fdmConfig;`numberOfSpotSteps`numberOfTimeSteps!(100;1000)];
    .config.mergeConfig[fdmConfig;`numberOfSpotSteps`numberOfTimeSteps!(200;2000)]);
gridBenchmark:.timing.benchmarkPricingGrid[trade;marketData;bsModel;configList];
.testutil.assertEqual[count gridBenchmark;2;"grid benchmark rows"];

-1 "PASS test_timing: pricingMs=",string timedPrice`elapsedMs;
