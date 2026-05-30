/ test_implied_vol_call.q - recover known call IV via bisection
\l core/init.q

trueVolatility:0.2;
marketPrice:.validation.blackScholesClosedForm[`call;100f;100f;1f;0.05;0f;trueVolatility];

ivResult:.iv.impliedVolatility[`call;marketPrice;100f;100f;1f;0.05;0f];

if[not ivResult[`status]~`OK; '"FAIL: IV solver status not OK"];
volError:abs ivResult[`impliedVolatility]-trueVolatility;
if[volError>1e-5; '"FAIL: implied vol error too large: ",string volError];
if[(abs ivResult`pricingError)>1e-6; '"FAIL: pricing error too large"];

-1 "PASS test_implied_vol_call: impliedVol=",string[ivResult`impliedVolatility],", trueVol=",string[trueVolatility],", error=",string volError;
