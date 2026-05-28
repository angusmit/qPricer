/ test_cache.q - pricing cache
\l ../../lib/init.q

/ Setup
trade1:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1f;`none;0Nf;0f);
mkt1:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
mdl:.model.createBlackScholesModel[];
cfg:.config.defaultPricingConfig[];

/ 1. Empty cache
emptyC:.cache.emptyCache[];
.testutil.assertTrue[0=.cache.cacheSize emptyC;"empty cache has 0 entries"];

/ 2. Same key for same inputs
keyA:.cache.createPricingKey[trade1;mkt1;cfg];
keyB:.cache.createPricingKey[trade1;mkt1;cfg];
.testutil.assertTrue[keyA~keyB;"same inputs produce same key"];

/ 3. Different key when spot changes
mkt2:@[mkt1;`spot;:;101f];
keyC:.cache.createPricingKey[trade1;mkt2;cfg];
.testutil.assertTrue[not keyA~keyC;"different spot produces different key"];

/ 4. Different key when vol changes
mkt3:@[mkt1;`volatility;:;0.25];
keyD:.cache.createPricingKey[trade1;mkt3;cfg];
.testutil.assertTrue[not keyA~keyD;"different vol produces different key"];

/ 5. Put then get
fakeResult:`unitPrice`testField!(42f;`test);
filledCache:.cache.put[emptyC;keyA;fakeResult];
retrieved:.cache.get[filledCache;keyA];
.testutil.assertTrue[retrieved[`unitPrice]=42f;"get returns stored result"];

/ 6. Contains
.testutil.assertTrue[.cache.contains[filledCache;keyA];"contains returns true after put"];
.testutil.assertTrue[not .cache.contains[filledCache;keyC];"contains returns false for missing"];

/ 7. getOrPrice matches normal pricing
normalResult:.engine.priceOption[trade1;mkt1;mdl;cfg];
cacheResult:.cache.getOrPrice[emptyC;trade1;mkt1;mdl;cfg];
.testutil.assertNear[cacheResult[`pricingResult][`unitPrice];normalResult`unitPrice;1e-12;"getOrPrice matches normal"];

/ 8. Local-vol key returns NOCACHE
lvMkt:.market.createLocalVolatilityMarketData[`AAPL;100f;0.05;0f;{[spotValue;timePoint] 0.2}];
lvKey:.cache.createPricingKey[trade1;lvMkt;cfg];
.testutil.assertTrue[lvKey~`NOCACHE;"local-vol key is NOCACHE"];

-1 "PASS test_cache";
