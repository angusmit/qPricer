/ test_local_vol_time_dependence.q - solver uses timePoint argument
\l core/init.q

callTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1f);

constantVolFunction:{[spotValue;timePoint] 0.2};
timeDependentFunction:{[spotValue;timePoint] 0.15+0.10*timePoint};

constantVolMkt:.market.createLocalVolatilityMarketData[`AAPL;100f;0.05;0f;constantVolFunction];
timeDependentMkt:.market.createLocalVolatilityMarketData[`AAPL;100f;0.05;0f;timeDependentFunction];

localVolModel:.model.createLocalVolatilityModel[];
config:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

constantPrice:(.engine.priceOption[callTrade;constantVolMkt;localVolModel;config])`unitPrice;
timeDependentPrice:(.engine.priceOption[callTrade;timeDependentMkt;localVolModel;config])`unitPrice;

/ 1. Both positive
if[constantPrice<=0f; '"FAIL: constant local-vol price not positive"];
if[timeDependentPrice<=0f; '"FAIL: time-dependent local-vol price not positive"];

/ 2. Prices differ (proves solver uses timePoint)
priceDifference:abs constantPrice-timeDependentPrice;
if[priceDifference<0.001; '"FAIL: time-dependent price did not differ from constant, solver may ignore timePoint"];

-1 "PASS test_local_vol_time_dependence: constant=",string[constantPrice],", timeDependent=",string[timeDependentPrice],", diff=",string priceDifference;
