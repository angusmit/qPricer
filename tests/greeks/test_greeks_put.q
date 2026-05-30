/ test_greeks_put.q — validate FDM Greeks for European put
\l core/init.q

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    2;`AAPL;`equityOption;`european;`put;100f;1f;1f);
marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(
    `AAPL;100f;0.05;0f;0.2);
model:.model.createBlackScholesModel[];
config:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

fdmG:.greeks.calculateGreeks[trade;marketData;model;config];
bsG:.validation.blackScholesGreeks[`put;100f;100f;1f;0.05;0f;0.2];

dErr:abs (first fdmG`delta)-bsG`delta;
gErr:abs (first fdmG`gamma)-bsG`gamma;
tErr:abs (first fdmG`theta)-bsG`theta;
vErr:abs (first fdmG`vega)-bsG`vega;
rErr:abs (first fdmG`rho)-bsG`rho;

allPass:(dErr<0.02) and (gErr<0.005) and (tErr<0.5) and (vErr<1.0) and rErr<1.0;

if[not allPass; '"FAIL test_greeks_put: d=",string[dErr]," g=",string[gErr]," t=",string[tErr]," v=",string[vErr]," r=",string rErr];
-1 "PASS test_greeks_put: delta=",string[first fdmG`delta],", gamma=",string[first fdmG`gamma],", theta=",string[first fdmG`theta],", vega=",string[first fdmG`vega],", rho=",string first fdmG`rho;
