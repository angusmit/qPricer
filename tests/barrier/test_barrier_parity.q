/ test_barrier_parity.q - knockIn + knockOut = vanilla
\l lib/init.q
mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
mdl:.model.createBlackScholesModel[];
cfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;300;4000;0f;300f;`linear;1b;1b);

/ Test up barrier call parity
vanillaCall:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`barrierType`barrierLevel`rebate!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1f;`none;0Nf;0f);
upOutCall:@[vanillaCall;`barrierType;:;`upAndOut]; upOutCall:@[upOutCall;`barrierLevel;:;130f];
upInCall:@[vanillaCall;`barrierType;:;`upAndIn]; upInCall:@[upInCall;`barrierLevel;:;130f];

vPrice:(.engine.priceOption[vanillaCall;mkt;mdl;cfg])`unitPrice;
koPrice:(.engine.priceOption[upOutCall;mkt;mdl;cfg])`unitPrice;
kiPrice:(.engine.priceOption[upInCall;mkt;mdl;cfg])`unitPrice;

parityError:abs (kiPrice+koPrice)-vPrice;
if[parityError>0.01; '"FAIL: up barrier parity error=",string parityError];
-1 "  PASS up barrier call parity: vanilla=",string[vPrice],", KO+KI=",string kiPrice+koPrice;

/ Test down barrier put parity
vanillaPut:@[vanillaCall;`optionType;:;`put];
dnOutPut:@[vanillaPut;`barrierType;:;`downAndOut]; dnOutPut:@[dnOutPut;`barrierLevel;:;70f];
dnInPut:@[vanillaPut;`barrierType;:;`downAndIn]; dnInPut:@[dnInPut;`barrierLevel;:;70f];

vpPrice:(.engine.priceOption[vanillaPut;mkt;mdl;cfg])`unitPrice;
koPrice2:(.engine.priceOption[dnOutPut;mkt;mdl;cfg])`unitPrice;
kiPrice2:(.engine.priceOption[dnInPut;mkt;mdl;cfg])`unitPrice;

parityError2:abs (kiPrice2+koPrice2)-vpPrice;
if[parityError2>0.01; '"FAIL: down barrier parity error=",string parityError2];
-1 "  PASS down barrier put parity: vanilla=",string[vpPrice],", KO+KI=",string kiPrice2+koPrice2;

-1 "PASS test_barrier_parity";
