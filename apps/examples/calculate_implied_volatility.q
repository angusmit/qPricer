/ calculate_implied_volatility.q - implied vol recovery example
/ Usage: q examples/calculate_implied_volatility.q

\l core/init.q
-1 "qFDM v",.qfdm.version," - Implied Volatility Example\n";

trueVol:0.2;
marketPrice:.validation.blackScholesClosedForm[`call;100f;100f;1f;0.05;0f;trueVol];
-1 "Market price from BS vol 20%: ",string marketPrice;

ivResult:.iv.impliedVolatility[`call;marketPrice;100f;100f;1f;0.05;0f];
-1 "Recovered implied volatility: ",string ivResult`impliedVolatility;
-1 "Pricing error:               ",string ivResult`pricingError;
-1 "Iterations:                  ",string ivResult`iterations;
-1 "";

/ Put example
putPrice:.validation.blackScholesClosedForm[`put;100f;100f;1f;0.05;0f;trueVol];
putIV:.iv.impliedVolatility[`put;putPrice;100f;100f;1f;0.05;0f];
-1 "Put market price:            ",string putPrice;
-1 "Put implied volatility:      ",string putIV`impliedVolatility;
