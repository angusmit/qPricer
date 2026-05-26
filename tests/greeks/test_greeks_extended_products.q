/ test_greeks_extended_products.q - Greeks for expanded product set
\l lib/init.q
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0.01;0.2);
mdl:.model.createBlackScholesModel[];
cfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

/ American call
amCall:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`american;`call;100f;1f;1f);
amCallGreeks:.greeks.calculateGreeks[amCall;mkt;mdl;cfg];
if[not (amCallGreeks[`delta]0)>0f; '"FAIL: American call delta should be positive"];
-1 "  PASS American call Greeks";

/ American put
amPut:@[amCall;`optionType;:;`put];
amPutGreeks:.greeks.calculateGreeks[amPut;mkt;mdl;cfg];
if[not (amPutGreeks[`delta]0)<0f; '"FAIL: American put delta should be negative"];
-1 "  PASS American put Greeks";

/ Barrier call (upAndOut)
barrierCall:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
    3;`AAPL;`equityOption;`european;`call;100f;1f;1f;`upAndOut;130f;0f);
barrierGreeks:.greeks.calculateGreeks[barrierCall;mkt;mdl;cfg];
if[not (barrierGreeks[`delta]0)>0f; '"FAIL: barrier call delta should be positive"];
-1 "  PASS Barrier call Greeks";

/ Knock-in call
knockInCall:@[barrierCall;`barrierType;:;`upAndIn];
kiGreeks:.greeks.calculateGreeks[knockInCall;mkt;mdl;cfg];
/ Knock-in delta might be any sign, just check it computes
if[null (kiGreeks[`delta]0); '"FAIL: knock-in delta should not be null"];
-1 "  PASS Knock-in call Greeks";

-1 "PASS test_greeks_extended_products";
