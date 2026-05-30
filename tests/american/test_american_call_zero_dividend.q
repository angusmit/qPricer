/ test_american_call_zero_dividend.q
\l core/init.q
amTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`american;`call;100f;1f;1f);
euTrade:@[amTrade;`exerciseStyle;:;`european];
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
mdl:.model.createBlackScholesModel[];
cfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);
amPrice:(.engine.priceOption[amTrade;mkt;mdl;cfg])`unitPrice;
euPrice:(.engine.priceOption[euTrade;mkt;mdl;cfg])`unitPrice;
priceDiff:abs amPrice-euPrice;
if[priceDiff>0.01; '"FAIL: American call (q=0) should approx equal European, diff=",string priceDiff];
-1 "PASS test_american_call_zero_dividend: american=",string[amPrice],", european=",string[euPrice],", diff=",string priceDiff;
