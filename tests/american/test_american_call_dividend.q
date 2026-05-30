/ test_american_call_dividend.q - American call with dividend >= European call
\l core/init.q
amTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`american;`call;100f;1f;1f);
euTrade:@[amTrade;`exerciseStyle;:;`european];
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0.03;0.2);
mdl:.model.createBlackScholesModel[];
cfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);
amPrice:(.engine.priceOption[amTrade;mkt;mdl;cfg])`unitPrice;
euPrice:(.engine.priceOption[euTrade;mkt;mdl;cfg])`unitPrice;
earlyExPremium:amPrice-euPrice;
if[earlyExPremium< -0.001; '"FAIL: American call with dividend should be >= European"];
-1 "PASS test_american_call_dividend: american=",string[amPrice],", european=",string[euPrice],", premium=",string earlyExPremium;
