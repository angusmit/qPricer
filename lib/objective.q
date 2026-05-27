/ objective.q - calibration objective functions (v0.21)

.objective.squaredError:{[modelPrice;marketPrice]
    diff:modelPrice-marketPrice;
    diff*diff
 };

.objective.absoluteError:{[modelPrice;marketPrice] abs modelPrice-marketPrice};

.objective.relativeError:{[modelPrice;marketPrice]
    if[marketPrice=0f; :0Nf];
    abs (modelPrice-marketPrice)%marketPrice
 };

.objective.rmse:{[modelPriceVector;marketPriceVector]
    diffs:modelPriceVector-marketPriceVector;
    sqrt avg diffs*diffs
 };

.objective.mae:{[modelPriceVector;marketPriceVector]
    avg abs modelPriceVector-marketPriceVector
 };

.objective.meanRelativeError:{[modelPriceVector;marketPriceVector]
    relErrors:{[mp;mk] $[mk=0f;0Nf;abs (mp-mk)%mk]}'[modelPriceVector;marketPriceVector];
    avg relErrors where not null relErrors
 };

.objective.weightedSse:{[modelPriceVector;marketPriceVector;weightVector]
    diffs:modelPriceVector-marketPriceVector;
    sum weightVector*diffs*diffs
 };

.objective.calibrationSummary:{[resultTable]
    / Use column access (not each) since resultTable may be a proper table
    modelPrices:resultTable`modelPrice;
    marketPrices:resultTable`marketPrice;
    absErrors:resultTable`absoluteError;
    `rowCount`rmse`mae`meanRelativeError`maxAbsoluteError`totalSse!(
        count resultTable;
        .objective.rmse[modelPrices;marketPrices];
        .objective.mae[modelPrices;marketPrices];
        .objective.meanRelativeError[modelPrices;marketPrices];
        max absErrors;
        sum {x*x} each modelPrices-marketPrices)
 };
