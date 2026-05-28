/ commoditySpread.q - spread option pricing (Kirk approximation) (v0.32)

/ params dict keys: fwd1, fwd2, strike, expiry, vol1, vol2, correlation, riskFreeRate
.commodity.spread.validateInputs:{[optType;params]
    if[params[`fwd1]<=0f; '"Non-positive forward 1"];
    if[params[`fwd2]<=0f; '"Non-positive forward 2"];
    if[params[`expiry]<=0f; '"Non-positive expiry"];
    if[params[`vol1]<=0f; '"Non-positive vol1"];
    if[params[`vol2]<=0f; '"Non-positive vol2"];
    corrVal:params`correlation;
    if[(corrVal< -1f) or corrVal>1f; '"Correlation out of range"];
    if[not optType in `call`put; '"Invalid option type"];
 };

.commodity.spread.spreadPayoff:{[fwd1;fwd2;strikePrice;optType]
    spreadVal:fwd1-fwd2;
    $[optType=`call;0f|spreadVal-strikePrice;0f|strikePrice-spreadVal]
 };

/ Kirk approximation for spread option
/ params dict keys: fwd1, fwd2, strike, expiry, vol1, vol2, correlation, riskFreeRate
.commodity.spread.kirkPrice:{[optType;params]
    .commodity.spread.validateInputs[optType;params];
    fwd1:params`fwd1; fwd2:params`fwd2;
    strikePrice:params`strike; expiryVal:params`expiry;
    vol1:params`vol1; vol2:params`vol2;
    corrVal:params`correlation; rateVal:params`riskFreeRate;
    adjStrike:fwd2+strikePrice;
    if[adjStrike<=0f; '"Non-positive adjusted strike"];
    f2Ratio:fwd2%adjStrike;
    vol2Adj:vol2*f2Ratio;
    twoCorr:2f*corrVal*vol1*vol2Adj;
    kirkVol:sqrt (vol1*vol1)+(vol2Adj*vol2Adj)-twoCorr;
    .commodity.black76.price[optType;fwd1;adjStrike;expiryVal;kirkVol;rateVal]
 };

/ Calendar spread option (long near, short far)
.commodity.spread.calendarSpreadPrice:{[optType;params]
    .commodity.spread.kirkPrice[optType;params]
 };
