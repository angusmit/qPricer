/ test_local_vol_flat_equivalence_put.q - flat local vol = flat BS for put
\l lib/init.q

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    2;`AAPL;`equityOption;`european;`put;100f;1f;1f);

flatMkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
localVolFunction:{[spotValue;timePoint] 0.2};
localVolMkt:.market.createLocalVolatilityMarketData[`AAPL;100f;0.05;0f;localVolFunction];

bsModel:.model.createBlackScholesModel[];
lvModel:.model.createLocalVolatilityModel[];

config:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

flatResult:.engine.priceOption[trade;flatMkt;bsModel;config];
localVolResult:.engine.priceOption[trade;localVolMkt;lvModel;config];
flatUnitPrice:flatResult`unitPrice;
localVolUnitPrice:localVolResult`unitPrice;
priceDiff:abs flatUnitPrice-localVolUnitPrice;
tolerance:0.01;

if[priceDiff>tolerance; '"FAIL test_local_vol_flat_equivalence_put: diff=",string priceDiff];
if[not localVolUnitPrice>0f; '"FAIL: local vol put price not positive"];

-1 "PASS test_local_vol_flat_equivalence_put: flat=",string[flatUnitPrice],", localVol=",string[localVolUnitPrice],", diff=",string priceDiff;
