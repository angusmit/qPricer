/ test_local_vol_skew_sanity.q - non-flat local vol changes prices sensibly
/ Uses a downside-only local volatility uplift: vol >= 0.2 everywhere,
/ with higher volatility below strike and flat volatility above strike.
\l core/init.q

callTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1f);
putTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    2;`AAPL;`equityOption;`european;`put;100f;1f;1f);

flatMarketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);

downsideOnlySkewFunction:{[spotValue;timePoint]
    extraVol:0f|0.001*(100f-spotValue);
    0.2+extraVol
 };

/ Verify function shape before pricing
lowSpotVolatility:downsideOnlySkewFunction[80f;0f];
atTheMoneyVolatility:downsideOnlySkewFunction[100f;0f];
highSpotVolatility:downsideOnlySkewFunction[120f;0f];

if[lowSpotVolatility<=atTheMoneyVolatility; '"FAIL: lower spot should have higher local volatility"];
if[highSpotVolatility<>atTheMoneyVolatility; '"FAIL: higher spot should equal flat volatility in downside-only skew"];
if[atTheMoneyVolatility<>0.2; '"FAIL: ATM local volatility should equal flat volatility"];

localVolMarketData:.market.createLocalVolatilityMarketData[`AAPL;100f;0.05;0f;downsideOnlySkewFunction];

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

/ 3. Vol >= flat everywhere, so both prices should increase
if[localVolCallPrice<flatCallPrice; '"FAIL: vol >= flat everywhere so call should not decrease"];
if[localVolPutPrice<flatPutPrice; '"FAIL: vol >= flat everywhere so put should not decrease"];

-1 "PASS test_local_vol_skew_sanity: flatCall=",string[flatCallPrice],", lvCall=",string[localVolCallPrice],", flatPut=",string[flatPutPrice],", lvPut=",string localVolPutPrice;
