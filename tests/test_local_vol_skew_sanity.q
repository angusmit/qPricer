/ test_local_vol_skew_sanity.q - non-flat local vol changes prices sensibly
\l lib/init.q

callTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1f);
putTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    2;`AAPL;`equityOption;`european;`put;100f;1f;1f);

flatMarketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);

downsideSkewFunction:{[spotValue;timePoint]
    baseVolatility:0.2;
    skewAdjustment:0.001*(100f-spotValue);
    candidateVolatility:baseVolatility+skewAdjustment;
    0.05|candidateVolatility&0.80
 };

localVolMarketData:.market.createLocalVolatilityMarketData[`AAPL;100f;0.05;0f;downsideSkewFunction];

blackScholesModel:.model.createBlackScholesModel[];
localVolModel:.model.createLocalVolatilityModel[];

config:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

flatCallPrice:(.engine.priceOption[callTrade;flatMarketData;blackScholesModel;config])`unitPrice;
localVolCallPrice:(.engine.priceOption[callTrade;localVolMarketData;localVolModel;config])`unitPrice;

flatPutPrice:(.engine.priceOption[putTrade;flatMarketData;blackScholesModel;config])`unitPrice;
localVolPutPrice:(.engine.priceOption[putTrade;localVolMarketData;localVolModel;config])`unitPrice;

/ 1. Positive prices
if[localVolCallPrice<=0f; '"FAIL: local-vol call price not positive"];
if[localVolPutPrice<=0f; '"FAIL: local-vol put price not positive"];

/ 2. Prices differ from flat
differenceTolerance:0.001;
callPriceDifference:abs localVolCallPrice-flatCallPrice;
putPriceDifference:abs localVolPutPrice-flatPutPrice;

if[callPriceDifference<differenceTolerance; '"FAIL: local-vol call price did not change versus flat vol"];
if[putPriceDifference<differenceTolerance; '"FAIL: local-vol put price did not change versus flat vol"];

-1 "PASS test_local_vol_skew_sanity: flatCall=",string[flatCallPrice],", lvCall=",string[localVolCallPrice],", flatPut=",string[flatPutPrice],", lvPut=",string localVolPutPrice;
