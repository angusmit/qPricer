/ strategy.q - generic strategy/backtest engine + path adapters + concrete strategies (v0.42)
/ -----------------------------------------------------------------------------------------
/ Layered design:
/   1. .strategy.*           registry-based generic engine; no strategy-specific code in .strategy.run
/   2. .strategy.path.*      data-source-agnostic path adapters (synthetic + Barchart)
/   3. .strategy.gammaScalp.* delta-hedged gamma scalping; self-registers at load
/ -----------------------------------------------------------------------------------------
/ Strategy interface contract (each strategy provides these four functions):
/   initFn[trade;firstStepRowDict;model;fdmConfig;stratCfg]        -> state dict
/   stepFn[state;marketStepRowDict;trade;model;fdmConfig;stratCfg] -> state dict
/   summaryFn[resultTable;stratCfg]                                -> summary dict
/   defaultCfgFn[]                                                 -> default stratCfg dict
/ -----------------------------------------------------------------------------------------
/ State convention: every state dict must carry the key rowEmit (a dict whose keys are the
/ per-step result-table columns). The driver folds stepFn over the path rows via scan and
/ then extracts the rowEmit dicts and builds the result table column-wise once.
/ -----------------------------------------------------------------------------------------
/ Path table standard schema (one row per step):
/   stepIndex (long), stepDate (date), spot (float), volatility (float),
/   riskFreeRate (float), dividendYield (float), marketPrice (float), status (symbol)
/ -----------------------------------------------------------------------------------------
/ Result table schema for gammaScalp extends the spec with two extra columns needed for
/ the gammaReconResidual cross-check and impliedVolAtEntry summary stat:
/   stepIndex, stepDate, spot, volatility, optionPrice, delta, hedgePosition, hedgeTrade,
/   txnCost, optionPnl, hedgePnl, financingPnl, thetaPnl, stepPnl, cumulativePnl,
/   theoreticalGammaPnl, status, message

/ ==================================================================
/ 1. Generic registry-based engine
/ ==================================================================

.strategy.registry:()!();

.strategy.register:{[strategyName;initFn;stepFn;summaryFn;defaultCfgFn]
    if[not -11h=type strategyName; '"strategy.register: strategyName must be a symbol atom"];
    funcs:`initFn`stepFn`summaryFn`defaultCfgFn!(initFn;stepFn;summaryFn;defaultCfgFn);
    .strategy.registry,:enlist[strategyName]!enlist funcs;
    strategyName
 };

.strategy.registeredStrategies:{[] key .strategy.registry};

.strategy.__lookup:{[strategyName]
    if[not strategyName in key .strategy.registry;
        '"strategy not registered: ",string strategyName];
    .strategy.registry strategyName
 };

.strategy.defaultConfig:{[strategyName]
    funcs:.strategy.__lookup strategyName;
    (funcs`defaultCfgFn)[]
 };

/ Column-wise build of a results table from a list of identically-shaped row dicts.
.strategy.__rowDictsToTable:{[rowDicts]
    if[0=count rowDicts; :()];
    colNames:key first rowDicts;
    flip colNames!{[colKey;rowsLocal] rowsLocal[;colKey]}[;rowDicts] each colNames
 };

/ Wrap stepFn with per-step try-catch: failure -> ERROR row, run continues.
.strategy.__safeStep:{[stepFnLocal;trade;model;fdmConfig;stratCfg;state;marketStep]
    result:.[stepFnLocal;(state;marketStep;trade;model;fdmConfig;stratCfg);{x}];
    if[10h=type result;
        errRowEmit:@[state`rowEmit;`status`message;:;(`ERROR;result)];
        errRowEmit:@[errRowEmit;`stepIndex;:;marketStep`stepIndex];
        errRowEmit:@[errRowEmit;`stepDate;:;marketStep`stepDate];
        errRowEmit:@[errRowEmit;`spot;:;marketStep`spot];
        :@[state;`rowEmit;:;errRowEmit]];
    result
 };

.strategy.run:{[strategyName;trade;pathTable;model;fdmConfig;stratCfg]
    funcs:.strategy.__lookup strategyName;
    .strategy.path.validate pathTable;
    pathRows:0!pathTable;
    firstStep:first pathRows;
    initialState:.[funcs`initFn;(trade;firstStep;model;fdmConfig;stratCfg);{x}];
    if[10h=type initialState; '"strategy.run init failed: ",initialState];
    remainingRows:1_pathRows;
    stepFnLocal:funcs`stepFn;
    safeFolder:.strategy.__safeStep[stepFnLocal;trade;model;fdmConfig;stratCfg;;];
    foldedStates:safeFolder\[initialState;remainingRows];
    states:enlist[initialState],foldedStates;
    rowDicts:states[;`rowEmit];
    .strategy.__rowDictsToTable rowDicts
 };

.strategy.summarize:{[strategyName;resultTable;stratCfg]
    funcs:.strategy.__lookup strategyName;
    (funcs`summaryFn)[resultTable;stratCfg]
 };

.strategy.runAndSummarize:{[strategyName;trade;pathTable;model;fdmConfig;stratCfg]
    resultTable:.strategy.run[strategyName;trade;pathTable;model;fdmConfig;stratCfg];
    summaryDict:.strategy.summarize[strategyName;resultTable;stratCfg];
    `result`summary!(resultTable;summaryDict)
 };

/ ==================================================================
/ 2. Path adapters
/ ==================================================================

.strategy.path.__schemaCols:`stepIndex`stepDate`spot`volatility`riskFreeRate`dividendYield`marketPrice`status;

.strategy.path.schema:{[]
    ([] stepIndex:0#0; stepDate:0#0Nd; spot:0#0Nf; volatility:0#0Nf; riskFreeRate:0#0Nf; dividendYield:0#0Nf; marketPrice:0#0Nf; status:0#`)
 };

.strategy.path.validate:{[pathTable]
    if[not 98h=type pathTable; '"strategy.path.validate: expected a table"];
    if[0=count pathTable; '"strategy.path.validate: empty path table"];
    missing:.strategy.path.__schemaCols where not .strategy.path.__schemaCols in cols pathTable;
    if[0<count missing; '"strategy.path.validate: missing columns ",", " sv string missing];
 };

/ Generate a deterministic GBM spot path. Total rows = steps; first row is spot0 at
/ stepIndex 0, remaining steps-1 rows are log-Brownian increments.
.strategy.path.fromSynthetic:{[pathCfg]
    required:`spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed;
    missing:required where not required in key pathCfg;
    if[0<count missing; '"strategy.path.fromSynthetic missing keys: ",", " sv string missing];
    if[pathCfg[`steps]<=1; '"strategy.path.fromSynthetic steps must be > 1"];
    if[pathCfg[`stepYears]<=0f; '"strategy.path.fromSynthetic stepYears must be positive"];
    if[pathCfg[`volatility]<0f; '"strategy.path.fromSynthetic volatility must be non-negative"];
    if[pathCfg[`spot0]<=0f; '"strategy.path.fromSynthetic spot0 must be positive"];
    stepsTotal:pathCfg`steps;
    spot0:`float$pathCfg`spot0;
    drift:`float$pathCfg`drift;
    vol:`float$pathCfg`volatility;
    dt:`float$pathCfg`stepYears;
    seed:pathCfg`seed;
    rfr:`float$pathCfg`riskFreeRate;
    divY:`float$pathCfg`dividendYield;
    normals:.montecarlo.__generateNormals[stepsTotal-1;seed];
    logIncrements:((drift-(0.5*vol*vol))*dt)+(vol*sqrt dt)*normals;
    logSpots:(log spot0)+sums logIncrements;
    spots:spot0,exp logSpots;
    startDate:$[`startDate in key pathCfg; pathCfg`startDate; 2024.01.01];
    stepDates:startDate+til stepsTotal;
    flip .strategy.path.__schemaCols!(
        til stepsTotal;
        stepDates;
        spots;
        stepsTotal#vol;
        stepsTotal#rfr;
        stepsTotal#divY;
        stepsTotal#0Nf;
        stepsTotal#`OK)
 };

/ Adapt an existing normalised Barchart contract series to the standard path schema.
/ No file IO; consumes the output of .parser.barchart.normalise only.
.strategy.path.fromBarchart:{[normalisedTable;contractId]
    if[not 98h=type normalisedTable; '"strategy.path.fromBarchart: expected normalised table"];
    if[not `contractId in cols normalisedTable;
        '"strategy.path.fromBarchart: normalisedTable missing contractId column"];
    contractRows:normalisedTable where (normalisedTable`contractId)=contractId;
    if[0=count contractRows; '"strategy.path.fromBarchart: no rows for contractId ",string contractId];
    sortedRows:`snapshotDate xasc contractRows;
    rowCount:count sortedRows;
    flip .strategy.path.__schemaCols!(
        til rowCount;
        sortedRows`snapshotDate;
        `float$sortedRows`spot;
        `float$sortedRows`impliedVolatility;
        rowCount#0.05;
        rowCount#0f;
        `float$sortedRows`marketPrice;
        sortedRows`status)
 };

/ ==================================================================
/ 3. Concrete strategy: delta-hedged gamma scalping
/ ==================================================================

.strategy.gammaScalp.defaultConfig:{[]
    `optionSide`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears`useModelDelta!(
        `long;
        `interval;
        1;
        0.05;
        0f;
        0f;
        1f%252f;
        1b)
 };

/ Column-key order for per-step row emit dicts; must match init and step output.
.strategy.gammaScalp.__rowEmitCols:`stepIndex`stepDate`spot`volatility`optionPrice`delta`hedgePosition`hedgeTrade`txnCost`optionPnl`hedgePnl`financingPnl`thetaPnl`stepPnl`cumulativePnl`theoreticalGammaPnl`status`message;

.strategy.gammaScalp.__validateConfig:{[stratCfg]
    requiredKeys:`optionSide`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears`useModelDelta;
    missing:requiredKeys where not requiredKeys in key stratCfg;
    if[0<count missing; '"gammaScalp config missing keys: ",", " sv string missing];
    if[not stratCfg[`optionSide] in `long`short;
        '"gammaScalp optionSide must be long or short"];
    if[not stratCfg[`rebalanceMode] in `interval`band;
        '"gammaScalp rebalanceMode must be interval or band"];
    if[(stratCfg[`rebalanceMode]=`interval)&(stratCfg[`rebalanceInterval]<=0);
        '"gammaScalp rebalanceInterval must be positive in interval mode"];
    if[(stratCfg[`rebalanceMode]=`band)&(stratCfg[`deltaBand]<0f);
        '"gammaScalp deltaBand must be non-negative in band mode"];
    if[stratCfg[`stepYears]<=0f; '"gammaScalp stepYears must be positive"];
    if[not stratCfg[`useModelDelta]; '"gammaScalp useModelDelta must be 1b (vendor delta path not implemented)"];
 };

.strategy.gammaScalp.init:{[trade;firstStep;model;fdmConfig;stratCfg]
    .strategy.gammaScalp.__validateConfig stratCfg;
    optionUnits:trade[`notional]*$[(stratCfg`optionSide)=`long; 1f; -1f];
    spot:firstStep`spot;
    vol:firstStep`volatility;
    rfr:firstStep`riskFreeRate;
    divY:firstStep`dividendYield;
    mktData:.market.createFlatMarketData[trade`underlying;spot;rfr;divY;vol];
    priceRes:.engine.priceOption[trade;mktData;model;fdmConfig];
    greeksRes:.greeks.calculateGreeks[trade;mktData;model;fdmConfig];
    optionPrice:priceRes`unitPrice;
    deltaVal:first greeksRes`delta;
    gammaVal:first greeksRes`gamma;
    thetaVal:first greeksRes`theta;
    hedgePosition:neg optionUnits*deltaVal;
    hedgeTrade:hedgePosition;
    txnCost:(abs hedgeTrade)*spot*stratCfg`txnCostRate;
    initialStepPnl:neg txnCost;
    rowEmit:.strategy.gammaScalp.__rowEmitCols!(
        firstStep`stepIndex;firstStep`stepDate;spot;vol;optionPrice;deltaVal;
        hedgePosition;hedgeTrade;txnCost;
        0f;0f;0f;0f;
        initialStepPnl;initialStepPnl;0f;
        `OK;"");
    cash:neg ((hedgePosition*spot)+txnCost);
    `optionUnits`hedgePosition`hedgedDelta`cash`numRebalances`prevSpot`prevOptionPrice`prevGamma`prevTheta`cumulativePnl`rowEmit!(
        optionUnits;hedgePosition;deltaVal;cash;1;spot;optionPrice;gammaVal;thetaVal;initialStepPnl;rowEmit)
 };

.strategy.gammaScalp.step:{[state;marketStep;trade;model;fdmConfig;stratCfg]
    optionUnits:state`optionUnits;
    prevHedgePos:state`hedgePosition;
    hedgedDeltaPrev:state`hedgedDelta;
    cashPrev:state`cash;
    numRebalancesPrev:state`numRebalances;
    prevSpot:state`prevSpot;
    prevOptionPrice:state`prevOptionPrice;
    prevGamma:state`prevGamma;
    prevTheta:state`prevTheta;
    cumulativePnlPrev:state`cumulativePnl;

    spot:marketStep`spot;
    vol:marketStep`volatility;
    rfr:marketStep`riskFreeRate;
    divY:marketStep`dividendYield;
    stepIndexVal:marketStep`stepIndex;
    stepDateVal:marketStep`stepDate;

    stepYears:stratCfg`stepYears;
    rebalanceInterval:stratCfg`rebalanceInterval;
    deltaBand:stratCfg`deltaBand;
    rebalanceMode:stratCfg`rebalanceMode;
    txnCostRate:stratCfg`txnCostRate;
    financingRate:stratCfg`financingRate;

    mktData:.market.createFlatMarketData[trade`underlying;spot;rfr;divY;vol];
    priceRes:.engine.priceOption[trade;mktData;model;fdmConfig];
    greeksRes:.greeks.calculateGreeks[trade;mktData;model;fdmConfig];
    optionPrice:priceRes`unitPrice;
    deltaVal:first greeksRes`delta;
    gammaVal:first greeksRes`gamma;
    thetaVal:first greeksRes`theta;

    optionPnl:optionUnits*(optionPrice-prevOptionPrice);
    hedgePnl:prevHedgePos*(spot-prevSpot);
    financingPnl:(financingRate*cashPrev)*stepYears;
    thetaPnl:(optionUnits*prevTheta)*stepYears;
    spotMove:spot-prevSpot;
    theoreticalGammaPnl:(0.5*prevGamma*optionUnits)*(spotMove*spotMove);

    intervalTrigger:0=stepIndexVal mod rebalanceInterval;
    bandTrigger:(abs deltaVal-hedgedDeltaPrev)>deltaBand;
    shouldRebalance:$[rebalanceMode=`interval; intervalTrigger; bandTrigger];

    newHedgePos:$[shouldRebalance; neg optionUnits*deltaVal; prevHedgePos];
    hedgeTrade:newHedgePos-prevHedgePos;
    txnCost:(abs hedgeTrade)*spot*txnCostRate;
    newHedgedDelta:$[shouldRebalance; deltaVal; hedgedDeltaPrev];
    rebalanceIncrement:$[shouldRebalance&0<>hedgeTrade; 1; 0];
    numRebalances:numRebalancesPrev+rebalanceIncrement;

    stepPnl:(optionPnl+hedgePnl+financingPnl)-txnCost;
    cumulativePnl:cumulativePnlPrev+stepPnl;
    newCash:((cashPrev+financingPnl)-hedgeTrade*spot)-txnCost;

    rowEmit:.strategy.gammaScalp.__rowEmitCols!(
        stepIndexVal;stepDateVal;spot;vol;optionPrice;deltaVal;
        newHedgePos;hedgeTrade;txnCost;
        optionPnl;hedgePnl;financingPnl;thetaPnl;
        stepPnl;cumulativePnl;theoreticalGammaPnl;
        `OK;"");
    state,`hedgePosition`hedgedDelta`cash`numRebalances`prevSpot`prevOptionPrice`prevGamma`prevTheta`cumulativePnl`rowEmit!(
        newHedgePos;newHedgedDelta;newCash;numRebalances;spot;optionPrice;gammaVal;thetaVal;cumulativePnl;rowEmit)
 };

.strategy.gammaScalp.summary:{[resultTable;stratCfg]
    emptyResult:`strategyName`steps`totalPnl`optionPnlTotal`hedgePnlTotal`txnCostTotal`financingTotal`numRebalances`realizedVol`impliedVolAtEntry`theoreticalGammaPnlTotal`thetaPnlTotal`gammaReconResidual`maxDrawdown`meanStepPnl`stepPnlVol`status`errorMessage!(
        `gammaScalp;0;0Nf;0Nf;0Nf;0Nf;0Nf;0;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;`ERROR;"empty result table");
    if[(0=count resultTable)|not 98h=type resultTable; :emptyResult];
    okRowsMask:(resultTable`status)=`OK;
    okRows:resultTable where okRowsMask;
    stepCount:count resultTable;
    if[0=count okRows; :@[emptyResult;(`steps;`errorMessage);:;(stepCount;"no OK rows")]];

    totalsDict:first 0!select
        totalPnl:sum stepPnl,
        optionPnlTotal:sum optionPnl,
        hedgePnlTotal:sum hedgePnl,
        txnCostTotal:sum txnCost,
        financingTotal:sum financingPnl,
        theoreticalGammaPnlTotal:sum theoreticalGammaPnl,
        thetaPnlTotal:sum thetaPnl,
        meanStepPnl:avg stepPnl,
        stepPnlVol:dev stepPnl
        from okRows;

    numRebalancesVal:sum 0<>okRows`hedgeTrade;

    stepYears:stratCfg`stepYears;
    spots:okRows`spot;
    logReturns:$[1<count spots; 1_(log spots)-prev log spots; ()];
    nonNullReturns:logReturns where not null logReturns;
    realizedVol:$[(0<count nonNullReturns)&stepYears>0f;
        (dev nonNullReturns)%sqrt stepYears;
        0Nf];

    impliedVolAtEntry:first okRows`volatility;

    cumPnlSeries:sums okRows`stepPnl;
    drawdownSeries:(maxs cumPnlSeries)-cumPnlSeries;
    maxDrawdownVal:max drawdownSeries;

    pnlExclCosts:(totalsDict[`totalPnl]+totalsDict`txnCostTotal)-totalsDict`financingTotal;
    gammaReconResidual:pnlExclCosts-(totalsDict[`theoreticalGammaPnlTotal]+totalsDict`thetaPnlTotal);

    totalsDict,`strategyName`steps`numRebalances`realizedVol`impliedVolAtEntry`gammaReconResidual`maxDrawdown`status`errorMessage!(
        `gammaScalp;stepCount;numRebalancesVal;realizedVol;impliedVolAtEntry;gammaReconResidual;maxDrawdownVal;`OK;"")
 };

/ Self-register at load.
.strategy.register[
    `gammaScalp;
    .strategy.gammaScalp.init;
    .strategy.gammaScalp.step;
    .strategy.gammaScalp.summary;
    .strategy.gammaScalp.defaultConfig];
