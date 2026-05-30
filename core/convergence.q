/ convergence.q - Monte Carlo convergence diagnostics (v0.19)

/ --- Public ---

.convergence.mcConvergenceTable:{[pricingFn;pathCountList;benchmarkPrice]
    resultList:();
    loopIdx:0;
    while[loopIdx<count pathCountList;
        pathCountValue:pathCountList loopIdx;
        startTime:.z.p;
        priceResult:.[pricingFn;enlist pathCountValue;{x}];
        elapsedMs:(`long$.z.p-startTime)%1000000;
        convergenceRow:.convergence.__buildRow[pathCountValue;priceResult;benchmarkPrice;elapsedMs];
        resultList:resultList,enlist convergenceRow;
        loopIdx+:1];
    resultList
 };

.convergence.__buildRow:{[pathCountValue;priceResult;benchmarkPrice;elapsedMs]
    if[10h=type priceResult;
        :`pathCount`unitPrice`standardError`lowerConfidence`upperConfidence`benchmarkPrice`absoluteError`relativeError`containsBenchmark`elapsedMs`status`errorMessage!(
            pathCountValue;0Nf;0Nf;0Nf;0Nf;benchmarkPrice;0Nf;0Nf;0b;elapsedMs;`ERROR;priceResult)];
    unitPriceVal:$[`unitPrice in key priceResult;priceResult`unitPrice;priceResult`price];
    seVal:priceResult`standardError;
    lowerCI:priceResult`lowerConfidence;
    upperCI:priceResult`upperConfidence;
    absError:abs unitPriceVal-benchmarkPrice;
    relError:$[benchmarkPrice<>0f;absError%abs benchmarkPrice;0Nf];
    containsBench:(benchmarkPrice>=lowerCI) and benchmarkPrice<=upperCI;
    `pathCount`unitPrice`standardError`lowerConfidence`upperConfidence`benchmarkPrice`absoluteError`relativeError`containsBenchmark`elapsedMs`status`errorMessage!(
        pathCountValue;unitPriceVal;seVal;lowerCI;upperCI;benchmarkPrice;absError;relError;containsBench;elapsedMs;`OK;"")
 };

.convergence.compareMcToBenchmark:{[mcResult;benchmarkPrice]
    mcPrice:$[`unitPrice in key mcResult;mcResult`unitPrice;mcResult`price];
    absError:abs mcPrice-benchmarkPrice;
    relError:$[benchmarkPrice<>0f;absError%abs benchmarkPrice;0Nf];
    `mcPrice`benchmarkPrice`absoluteError`relativeError!(mcPrice;benchmarkPrice;absError;relError)
 };

.convergence.confidenceIntervalContains:{[mcResult;benchmarkPrice]
    lowerCI:mcResult`lowerConfidence;
    upperCI:mcResult`upperConfidence;
    (benchmarkPrice>=lowerCI) and benchmarkPrice<=upperCI
 };

.convergence.summariseConvergence:{[convergenceTable]
    okRows:convergenceTable where ({x[`status]~`OK} each convergenceTable);
    if[0=count okRows; :`minPathCount`maxPathCount`firstError`lastError`errorImproved`firstStandardError`lastStandardError`standardErrorImproved!(0N;0N;0Nf;0Nf;0b;0Nf;0Nf;0b)];
    pathCounts:{x`pathCount} each okRows;
    errors:{x`absoluteError} each okRows;
    stdErrors:{x`standardError} each okRows;
    firstError:errors 0;
    lastError:last errors;
    firstSE:stdErrors 0;
    lastSE:last stdErrors;
    `minPathCount`maxPathCount`firstError`lastError`errorImproved`firstStandardError`lastStandardError`standardErrorImproved!(
        min pathCounts;max pathCounts;firstError;lastError;lastError<=firstError;firstSE;lastSE;lastSE<=firstSE)
 };
