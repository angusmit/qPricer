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
/ 3. Shared delta-hedge accounting helpers
/ ==================================================================
/ Pure functions used by both gammaScalp and shortVariance. They take a generic
/ position view (signed positionValue, signed positionDelta) so the strategy
/ chooses whether the position is long/short, single-leg or multi-leg.
/ -----------------------------------------------------------------------------
/ __hedgeInit input dict (hedgeInputs) keys:
/   spot, positionDelta, txnCostRate
/ __hedgeInit output dict keys:
/   hedgePosition, hedgeTrade, txnCost, cashAdj
/ (cashAdj is the cash flow from setting up the hedge alone, not including any
/  premium received or paid for the underlying position itself.)
/ -----------------------------------------------------------------------------
/ __hedgeStep state (hedgeState) keys:
/   cash, hedgePosition, hedgedDelta, numRebalances, prevSpot, prevPositionValue
/ __hedgeStep input (stepInputs) keys:
/   spot, positionValue, positionDelta, stepIndex, stepYears,
/   txnCostRate, financingRate, rebalanceMode, rebalanceInterval, deltaBand
/ __hedgeStep output keys (updated state + per-step accounting):
/   cash, hedgePosition, hedgedDelta, numRebalances, prevSpot, prevPositionValue,
/   hedgeTrade, txnCost, financingPnl, hedgePnl, positionPnl, stepPnl

.strategy.__hedgeInit:{[hedgeInputs]
    spot:hedgeInputs`spot;
    positionDelta:hedgeInputs`positionDelta;
    txnCostRate:hedgeInputs`txnCostRate;
    hedgePosition:neg positionDelta;
    hedgeTrade:hedgePosition;
    txnCost:(abs hedgeTrade)*spot*txnCostRate;
    cashAdj:neg ((hedgePosition*spot)+txnCost);
    `hedgePosition`hedgeTrade`txnCost`cashAdj!(hedgePosition;hedgeTrade;txnCost;cashAdj)
 };

/ Total portfolio value from atomic components. Generic, strategy-agnostic.
/ Independent of any flow attribution; used to verify accounting tests.
.strategy.__portfolioValue:{[cashVal;legMarkSum;hedgePosition;spot]
    cashVal + legMarkSum + hedgePosition * spot
 };

/ Multi-underlying portfolio value: hedgePositions and spots are aligned vectors.
/ Additive helper; .strategy.__portfolioValue (single-asset) remains byte-identical.
.strategy.__portfolioValueMulti:{[cashVal;legMarkSum;hedgePositions;spots]
    cashVal + legMarkSum + sum hedgePositions * spots
 };

.strategy.__hedgeStep:{[hedgeState;stepInputs]
    spot:stepInputs`spot;
    positionValue:stepInputs`positionValue;
    positionDelta:stepInputs`positionDelta;
    stepIndexVal:stepInputs`stepIndex;
    stepYears:stepInputs`stepYears;
    txnCostRate:stepInputs`txnCostRate;
    financingRate:stepInputs`financingRate;
    rebalanceMode:stepInputs`rebalanceMode;
    rebalanceInterval:stepInputs`rebalanceInterval;
    deltaBand:stepInputs`deltaBand;

    prevSpot:hedgeState`prevSpot;
    prevPositionValue:hedgeState`prevPositionValue;
    prevHedgePos:hedgeState`hedgePosition;
    prevHedgedDelta:hedgeState`hedgedDelta;
    cashPrev:hedgeState`cash;
    numRebalancesPrev:hedgeState`numRebalances;

    positionPnl:positionValue-prevPositionValue;
    hedgePnl:prevHedgePos*(spot-prevSpot);
    financingPnl:(financingRate*cashPrev)*stepYears;

    intervalTrigger:0=stepIndexVal mod rebalanceInterval;
    bandTrigger:(abs positionDelta-prevHedgedDelta)>deltaBand;
    shouldRebalance:$[rebalanceMode=`interval; intervalTrigger; bandTrigger];

    newHedgePos:$[shouldRebalance; neg positionDelta; prevHedgePos];
    hedgeTrade:newHedgePos-prevHedgePos;
    txnCost:(abs hedgeTrade)*spot*txnCostRate;
    newHedgedDelta:$[shouldRebalance; positionDelta; prevHedgedDelta];
    rebalanceIncrement:$[shouldRebalance&0<>hedgeTrade; 1; 0];
    numRebalances:numRebalancesPrev+rebalanceIncrement;

    stepPnl:(positionPnl+hedgePnl+financingPnl)-txnCost;
    newCash:((cashPrev+financingPnl)-hedgeTrade*spot)-txnCost;

    `cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue`hedgeTrade`txnCost`financingPnl`hedgePnl`positionPnl`stepPnl!(
        newCash;newHedgePos;newHedgedDelta;numRebalances;spot;positionValue;hedgeTrade;txnCost;financingPnl;hedgePnl;positionPnl;stepPnl)
 };

/ ==================================================================
/ 4. Concrete strategy: delta-hedged gamma scalping
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
    positionValue:optionUnits*optionPrice;
    positionDelta:optionUnits*deltaVal;
    positionGamma:optionUnits*gammaVal;
    positionTheta:optionUnits*thetaVal;
    hedgeInit:.strategy.__hedgeInit `spot`positionDelta`txnCostRate!(spot;positionDelta;stratCfg`txnCostRate);
    initialStepPnl:neg hedgeInit`txnCost;
    rowEmit:.strategy.gammaScalp.__rowEmitCols!(
        firstStep`stepIndex;firstStep`stepDate;spot;vol;optionPrice;deltaVal;
        hedgeInit`hedgePosition;hedgeInit`hedgeTrade;hedgeInit`txnCost;
        0f;0f;0f;0f;
        initialStepPnl;initialStepPnl;0f;
        `OK;"");
    `optionUnits`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue`prevPositionGamma`prevPositionTheta`hedgeTrade`txnCost`financingPnl`hedgePnl`positionPnl`stepPnl`cumulativePnl`rowEmit!(
        optionUnits;hedgeInit`cashAdj;hedgeInit`hedgePosition;positionDelta;1;spot;positionValue;positionGamma;positionTheta;
        hedgeInit`hedgeTrade;hedgeInit`txnCost;0f;0f;0f;initialStepPnl;initialStepPnl;rowEmit)
 };

.strategy.gammaScalp.step:{[state;marketStep;trade;model;fdmConfig;stratCfg]
    spot:marketStep`spot;
    vol:marketStep`volatility;
    mktData:.market.createFlatMarketData[trade`underlying;spot;marketStep`riskFreeRate;marketStep`dividendYield;vol];
    priceRes:.engine.priceOption[trade;mktData;model;fdmConfig];
    greeksRes:.greeks.calculateGreeks[trade;mktData;model;fdmConfig];
    optionPrice:priceRes`unitPrice;
    deltaVal:first greeksRes`delta;
    gammaVal:first greeksRes`gamma;
    thetaVal:first greeksRes`theta;
    optionUnits:state`optionUnits;
    positionValue:optionUnits*optionPrice;
    positionDelta:optionUnits*deltaVal;
    positionGamma:optionUnits*gammaVal;
    positionTheta:optionUnits*thetaVal;

    hedgeState:`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue#state;
    stepInputs:`spot`positionValue`positionDelta`stepIndex`stepYears`txnCostRate`financingRate`rebalanceMode`rebalanceInterval`deltaBand!(
        spot;positionValue;positionDelta;marketStep`stepIndex;stratCfg`stepYears;stratCfg`txnCostRate;stratCfg`financingRate;stratCfg`rebalanceMode;stratCfg`rebalanceInterval;stratCfg`deltaBand);
    newHedge:.strategy.__hedgeStep[hedgeState;stepInputs];

    spotMove:spot-state`prevSpot;
    theoreticalGammaPnl:(0.5*state`prevPositionGamma)*spotMove*spotMove;
    thetaPnl:(state`prevPositionTheta)*stratCfg`stepYears;
    cumulativePnl:(state`cumulativePnl)+newHedge`stepPnl;

    rowEmit:.strategy.gammaScalp.__rowEmitCols!(
        marketStep`stepIndex;marketStep`stepDate;spot;vol;optionPrice;deltaVal;
        newHedge`hedgePosition;newHedge`hedgeTrade;newHedge`txnCost;
        newHedge`positionPnl;newHedge`hedgePnl;newHedge`financingPnl;thetaPnl;
        newHedge`stepPnl;cumulativePnl;theoreticalGammaPnl;
        `OK;"");
    state,newHedge,`prevPositionGamma`prevPositionTheta`cumulativePnl`rowEmit!(
        positionGamma;positionTheta;cumulativePnl;rowEmit)
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

/ ==================================================================
/ 5. Concrete strategy: short variance (variance risk premium)
/ ==================================================================
/ Sells a delta-hedged ATM straddle (one call + one put at trade.strike, same expiry)
/ when implied vol is rich vs a config-supplied realized-vol forecast. The straddle
/ is held flat to expiry; per-step hedging shares the .strategy.__hedgeStep helper
/ with gammaScalp. Premium is collected up front and added to cash; per-step PnL
/ attribution is option-MTM + hedge MTM + financing - txn cost, and is reconciled
/ against 0.5 * netGamma * (dSpot)^2 + netTheta * dt (both negative for short option).
/ -----------------------------------------------------------------------------
/ Assumption: the trade dict carries the SHORT straddle as a single trade. We build
/ the two legs internally by setting optionType=call/put on the same strike, expiry,
/ underlying, notional. The hedge accounts the SHORT position so positionValue is
/ negative; positionDelta is the net signed delta; positionGamma/Theta are negative
/ for a vanilla short straddle.
/ -----------------------------------------------------------------------------
/ Result-table schema extends the spec with `volatility, `thetaPnl, `premiumCollected
/ columns (per-step) so the summary can compute impliedVolAtEntry, varianceRiskPremium,
/ and gammaReconResidual from resultTable alone (summary receives no state).
/ Flat-path semantics: if the entry gate is closed (impliedVol <= forecastVol+entryMargin),
/ every step emits status=`flat with zero P&L and premiumCollected 0.

.strategy.shortVariance.__rowEmitCols:`stepIndex`stepDate`spot`volatility`callPrice`putPrice`positionValue`netDelta`hedgePosition`hedgeTrade`txnCost`positionPnl`hedgePnl`financingPnl`thetaPnl`stepPnl`cumulativePnl`theoreticalGammaPnl`premiumCollected`status`message;

.strategy.shortVariance.defaultConfig:{[]
    `forecastVol`entryMargin`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears!(
        0.20;
        0.02;
        `interval;
        1;
        0.05;
        0f;
        0f;
        1f%252f)
 };

.strategy.shortVariance.__validateConfig:{[stratCfg]
    requiredKeys:`forecastVol`entryMargin`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears;
    missing:requiredKeys where not requiredKeys in key stratCfg;
    if[0<count missing; '"shortVariance config missing keys: ",", " sv string missing];
    if[(stratCfg`forecastVol)<0f; '"shortVariance forecastVol must be non-negative"];
    if[(stratCfg`entryMargin)<0f; '"shortVariance entryMargin must be non-negative"];
    if[not stratCfg[`rebalanceMode] in `interval`band;
        '"shortVariance rebalanceMode must be interval or band"];
    if[(stratCfg[`rebalanceMode]=`interval)&(stratCfg[`rebalanceInterval]<=0);
        '"shortVariance rebalanceInterval must be positive in interval mode"];
    if[(stratCfg[`rebalanceMode]=`band)&(stratCfg[`deltaBand]<0f);
        '"shortVariance deltaBand must be non-negative in band mode"];
    if[(stratCfg`stepYears)<=0f; '"shortVariance stepYears must be positive"];
 };

/ Build the two straddle legs from the input trade.
.strategy.shortVariance.__buildLegs:{[trade]
    callTrade:@[trade;(`tradeId;`optionType);:;(`$(string trade`tradeId),"_C";`call)];
    putTrade:@[trade;(`tradeId;`optionType);:;(`$(string trade`tradeId),"_P";`put)];
    `callTrade`putTrade!(callTrade;putTrade)
 };

/ Build a flat-path row emit (gate closed) with zero PnL fields.
.strategy.shortVariance.__flatRow:{[marketStep]
    .strategy.shortVariance.__rowEmitCols!(
        marketStep`stepIndex;marketStep`stepDate;marketStep`spot;marketStep`volatility;
        0Nf;0Nf;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;
        `flat;"Entry gate closed")
 };

.strategy.shortVariance.init:{[trade;firstStep;model;fdmConfig;stratCfg]
    .strategy.shortVariance.__validateConfig stratCfg;
    notional:trade`notional;
    spot:firstStep`spot;
    vol:firstStep`volatility;
    rfr:firstStep`riskFreeRate;
    divY:firstStep`dividendYield;
    forecastVol:stratCfg`forecastVol;
    entryMargin:stratCfg`entryMargin;
    gateOpen:vol>forecastVol+entryMargin;
    legs:.strategy.shortVariance.__buildLegs trade;
    callTrade:legs`callTrade;
    putTrade:legs`putTrade;
    if[not gateOpen;
        flatRow:.strategy.shortVariance.__flatRow firstStep;
        flatRow:@[flatRow;`message;:;"Gate closed: vol ",(string vol)," <= forecast+margin ",string forecastVol+entryMargin];
        :`gateOpen`notional`callTrade`putTrade`impliedVolAtEntry`premiumCollected`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue`prevPositionGamma`prevPositionTheta`hedgeTrade`txnCost`financingPnl`hedgePnl`positionPnl`stepPnl`cumulativePnl`rowEmit!(
            0b;notional;callTrade;putTrade;vol;0f;0f;0f;0f;0;spot;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;flatRow)];
    mktData:.market.createFlatMarketData[trade`underlying;spot;rfr;divY;vol];
    callPrice:(.engine.priceOption[callTrade;mktData;model;fdmConfig])`unitPrice;
    putPrice:(.engine.priceOption[putTrade;mktData;model;fdmConfig])`unitPrice;
    callGreeks:.greeks.calculateGreeks[callTrade;mktData;model;fdmConfig];
    putGreeks:.greeks.calculateGreeks[putTrade;mktData;model;fdmConfig];
    callDelta:first callGreeks`delta;
    putDelta:first putGreeks`delta;
    callGamma:first callGreeks`gamma;
    putGamma:first putGreeks`gamma;
    callTheta:first callGreeks`theta;
    putTheta:first putGreeks`theta;
    positionValue:neg notional*callPrice+putPrice;
    positionDelta:neg notional*callDelta+putDelta;
    positionGamma:neg notional*callGamma+putGamma;
    positionTheta:neg notional*callTheta+putTheta;
    premiumCollected:notional*callPrice+putPrice;
    hedgeInit:.strategy.__hedgeInit `spot`positionDelta`txnCostRate!(spot;positionDelta;stratCfg`txnCostRate);
    initialCash:premiumCollected+hedgeInit`cashAdj;
    initialStepPnl:neg hedgeInit`txnCost;
    rowEmit:.strategy.shortVariance.__rowEmitCols!(
        firstStep`stepIndex;firstStep`stepDate;spot;vol;callPrice;putPrice;positionValue;positionDelta;
        hedgeInit`hedgePosition;hedgeInit`hedgeTrade;hedgeInit`txnCost;
        0f;0f;0f;0f;
        initialStepPnl;initialStepPnl;0f;
        premiumCollected;
        `OK;"");
    `gateOpen`notional`callTrade`putTrade`impliedVolAtEntry`premiumCollected`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue`prevPositionGamma`prevPositionTheta`hedgeTrade`txnCost`financingPnl`hedgePnl`positionPnl`stepPnl`cumulativePnl`rowEmit!(
        1b;notional;callTrade;putTrade;vol;premiumCollected;initialCash;hedgeInit`hedgePosition;positionDelta;1;spot;positionValue;positionGamma;positionTheta;
        hedgeInit`hedgeTrade;hedgeInit`txnCost;0f;0f;0f;initialStepPnl;initialStepPnl;rowEmit)
 };

.strategy.shortVariance.step:{[state;marketStep;trade;model;fdmConfig;stratCfg]
    if[not state`gateOpen;
        flatRow:.strategy.shortVariance.__flatRow marketStep;
        :@[state;`rowEmit;:;flatRow]];
    notional:state`notional;
    callTrade:state`callTrade;
    putTrade:state`putTrade;
    spot:marketStep`spot;
    vol:marketStep`volatility;
    mktData:.market.createFlatMarketData[trade`underlying;spot;marketStep`riskFreeRate;marketStep`dividendYield;vol];
    callPrice:(.engine.priceOption[callTrade;mktData;model;fdmConfig])`unitPrice;
    putPrice:(.engine.priceOption[putTrade;mktData;model;fdmConfig])`unitPrice;
    callGreeks:.greeks.calculateGreeks[callTrade;mktData;model;fdmConfig];
    putGreeks:.greeks.calculateGreeks[putTrade;mktData;model;fdmConfig];
    callDelta:first callGreeks`delta;
    putDelta:first putGreeks`delta;
    callGamma:first callGreeks`gamma;
    putGamma:first putGreeks`gamma;
    callTheta:first callGreeks`theta;
    putTheta:first putGreeks`theta;
    positionValue:neg notional*callPrice+putPrice;
    positionDelta:neg notional*callDelta+putDelta;
    positionGamma:neg notional*callGamma+putGamma;
    positionTheta:neg notional*callTheta+putTheta;
    hedgeState:`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue#state;
    stepInputs:`spot`positionValue`positionDelta`stepIndex`stepYears`txnCostRate`financingRate`rebalanceMode`rebalanceInterval`deltaBand!(
        spot;positionValue;positionDelta;marketStep`stepIndex;stratCfg`stepYears;stratCfg`txnCostRate;stratCfg`financingRate;stratCfg`rebalanceMode;stratCfg`rebalanceInterval;stratCfg`deltaBand);
    newHedge:.strategy.__hedgeStep[hedgeState;stepInputs];
    spotMove:spot-state`prevSpot;
    theoreticalGammaPnl:(0.5*state`prevPositionGamma)*spotMove*spotMove;
    thetaPnl:(state`prevPositionTheta)*stratCfg`stepYears;
    cumulativePnl:(state`cumulativePnl)+newHedge`stepPnl;
    rowEmit:.strategy.shortVariance.__rowEmitCols!(
        marketStep`stepIndex;marketStep`stepDate;spot;vol;callPrice;putPrice;positionValue;positionDelta;
        newHedge`hedgePosition;newHedge`hedgeTrade;newHedge`txnCost;
        newHedge`positionPnl;newHedge`hedgePnl;newHedge`financingPnl;thetaPnl;
        newHedge`stepPnl;cumulativePnl;theoreticalGammaPnl;
        state`premiumCollected;
        `OK;"");
    state,newHedge,`prevPositionGamma`prevPositionTheta`cumulativePnl`rowEmit!(
        positionGamma;positionTheta;cumulativePnl;rowEmit)
 };

.strategy.shortVariance.summary:{[resultTable;stratCfg]
    emptyResult:`strategyName`gateOpen`steps`premiumCollected`totalPnl`positionPnlTotal`hedgePnlTotal`txnCostTotal`financingTotal`numRebalances`impliedVolAtEntry`realizedVol`varianceRiskPremium`theoreticalGammaPnlTotal`thetaPnlTotal`gammaReconResidual`maxDrawdown`meanStepPnl`stepPnlVol`status`errorMessage!(
        `shortVariance;0b;0;0f;0Nf;0Nf;0Nf;0Nf;0Nf;0;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;`ERROR;"empty result table");
    if[(0=count resultTable)|not 98h=type resultTable; :emptyResult];
    statusCol:resultTable`status;
    stepCount:count resultTable;
    if[all statusCol=`flat;
        :`strategyName`gateOpen`steps`premiumCollected`totalPnl`positionPnlTotal`hedgePnlTotal`txnCostTotal`financingTotal`numRebalances`impliedVolAtEntry`realizedVol`varianceRiskPremium`theoreticalGammaPnlTotal`thetaPnlTotal`gammaReconResidual`maxDrawdown`meanStepPnl`stepPnlVol`status`errorMessage!(
            `shortVariance;0b;stepCount;0f;0f;0f;0f;0f;0f;0;first resultTable`volatility;0Nf;0Nf;0f;0f;0f;0f;0f;0f;`flat;"Entry gate closed")];
    okRowsMask:statusCol=`OK;
    okRows:resultTable where okRowsMask;
    if[0=count okRows; :@[emptyResult;`steps;:;stepCount]];

    totalsDict:first 0!select
        totalPnl:sum stepPnl,
        positionPnlTotal:sum positionPnl,
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
    varianceRiskPremium:impliedVolAtEntry-realizedVol;
    cumPnlSeries:sums okRows`stepPnl;
    drawdownSeries:(maxs cumPnlSeries)-cumPnlSeries;
    maxDrawdownVal:max drawdownSeries;
    pnlExclCosts:(totalsDict[`totalPnl]+totalsDict`txnCostTotal)-totalsDict`financingTotal;
    gammaReconResidual:pnlExclCosts-(totalsDict[`theoreticalGammaPnlTotal]+totalsDict`thetaPnlTotal);
    premiumCollected:first okRows`premiumCollected;

    totalsDict,`strategyName`gateOpen`steps`premiumCollected`numRebalances`impliedVolAtEntry`realizedVol`varianceRiskPremium`gammaReconResidual`maxDrawdown`status`errorMessage!(
        `shortVariance;1b;stepCount;premiumCollected;numRebalancesVal;impliedVolAtEntry;realizedVol;varianceRiskPremium;gammaReconResidual;maxDrawdownVal;`OK;"")
 };

.strategy.register[
    `shortVariance;
    .strategy.shortVariance.init;
    .strategy.shortVariance.step;
    .strategy.shortVariance.summary;
    .strategy.shortVariance.defaultConfig];

/ ==================================================================
/ 6. Concrete strategy: calendar roll (mutable leg lifecycle)
/ ==================================================================
/ Holds a 2-leg calendar spread (same strike, same option type, two expiries) and
/ rolls expiring legs to fresh tenors. First strategy whose position mutates mid-path.
/ -----------------------------------------------------------------------------
/ Accounting uses the portfolio-value identity, not __hedgeStep's positionPnl:
/   PV_t = cash + sum(legSide * legUnits * legMark) + hedgePosition * spot
/   stepPnl = PV_t - PV_{t-1}
/         == positionPnl (surviving-leg mark changes)
/          + rollPnl     (expiring-leg mark change from prev mark to settlement intrinsic)
/          + hedgePnl
/          + financingPnl
/          - txnCost     (sum of hedge txn cost + roll txn cost on closed + opened legs)
/ Cash flows: settlement intrinsic IN, new-leg entry mark OUT, roll txnCost OUT.
/ Mark-to-cash conversions at roll are PV-neutral; only the mark changes show as PnL.
/ -----------------------------------------------------------------------------
/ __hedgeStep is reused for the hedge leg only. Its returned positionPnl/stepPnl are
/ discarded; we use its hedge accounting (hedgeTrade/txnCost/hedgePnl/financingPnl/cash
/ update). Helper does not mutate the leg book.
/ -----------------------------------------------------------------------------
/ Expiry edge case: a leg with remainingTime <= rollThresholdYears is settled at
/ intrinsic via payoff (not priced via the FDM grid which is invalid at expiry).
/ Intrinsic delta is 1 (or -1 for puts) when ITM, 0 OTM, 0 ATM; intrinsic gamma/theta = 0.
/ -----------------------------------------------------------------------------
/ Result schema extends spec with thetaPnl column so non-roll-step reconciliation
/ can be computed from resultTable alone.
/ Result columns:
/   stepIndex, stepDate, spot, currentTimeYears, activeLegCount, frontRemainingTime,
/   backRemainingTime, positionValue, netDelta, hedgePosition, hedgeTrade, txnCost,
/   positionPnl, rollPnl, hedgePnl, financingPnl, thetaPnl, rollEvents, stepPnl,
/   cumulativePnl, theoreticalGammaPnl, status, message

.strategy.calendarRoll.__rowEmitCols:`stepIndex`stepDate`spot`currentTimeYears`activeLegCount`frontRemainingTime`backRemainingTime`positionValue`netDelta`hedgePosition`hedgeTrade`txnCost`positionPnl`rollPnl`hedgePnl`financingPnl`thetaPnl`rollEvents`stepPnl`cumulativePnl`theoreticalGammaPnl`status`message;

.strategy.calendarRoll.__emptyLegTable:{[]
    ([] legId:0#0; legRole:0#`; optionType:0#`; strike:0#0f; side:0#0f; units:0#0f;
        expiryTimeYears:0#0f; entryTimeYears:0#0f; entrySpot:0#0f; entryPrice:0#0f;
        currentMark:0#0f; currentDelta:0#0f; currentGamma:0#0f; currentTheta:0#0f)
 };

.strategy.calendarRoll.defaultConfig:{[]
    `spreadType`optionType`frontTenorYears`backTenorYears`rollThresholdYears`rollBackLeg`hedgeDelta`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears!(
        `longCalendar;
        `call;
        0.05;
        0.20;
        0.005;
        1b;
        1b;
        `interval;
        1;
        0.05;
        0f;
        0f;
        1f%252f)
 };

.strategy.calendarRoll.__validateConfig:{[stratCfg]
    requiredKeys:`spreadType`optionType`frontTenorYears`backTenorYears`rollThresholdYears`rollBackLeg`hedgeDelta`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears;
    missing:requiredKeys where not requiredKeys in key stratCfg;
    if[0<count missing; '"calendarRoll config missing keys: ",", " sv string missing];
    if[not stratCfg[`spreadType] in `longCalendar`shortCalendar;
        '"calendarRoll spreadType must be longCalendar or shortCalendar"];
    if[not stratCfg[`optionType] in `call`put;
        '"calendarRoll optionType must be call or put"];
    if[(stratCfg`frontTenorYears)<=0f; '"calendarRoll frontTenorYears must be positive"];
    if[(stratCfg`backTenorYears)<=stratCfg`frontTenorYears;
        '"calendarRoll backTenorYears must exceed frontTenorYears"];
    if[(stratCfg`rollThresholdYears)<0f; '"calendarRoll rollThresholdYears must be non-negative"];
    if[(stratCfg[`rebalanceMode]=`interval)&(stratCfg[`rebalanceInterval]<=0);
        '"calendarRoll rebalanceInterval must be positive in interval mode"];
    if[(stratCfg`stepYears)<=0f; '"calendarRoll stepYears must be positive"];
 };

/ Price (or settle at intrinsic) one leg at the given currentTimeYears. Returns a dict
/ with unitPrice, unitDelta, unitGamma, unitTheta, remainingTime, isExpiring.
.strategy.calendarRoll.__priceLeg:{[legRow;currentTimeYears;contextDict;model;fdmConfig]
    rt:legRow[`expiryTimeYears]-currentTimeYears;
    spot:contextDict`spot;
    rollThr:contextDict`rollThresholdYears;
    isExp:rt<=rollThr;
    if[rt<=1e-10;
        intrinsic:$[legRow[`optionType]=`call;
            0f|spot-legRow`strike;
            0f|(legRow`strike)-spot];
        intrinsicDelta:$[legRow[`optionType]=`call;
            $[spot>legRow`strike;1f;0f];
            $[spot<legRow`strike;-1f;0f]];
        :`unitPrice`unitDelta`unitGamma`unitTheta`remainingTime`isExpiring!(
            intrinsic;intrinsicDelta;0f;0f;rt;1b)];
    legTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
        `crleg;contextDict`underlying;`equityOption;`european;legRow`optionType;legRow`strike;rt;1f);
    mktData:.market.createFlatMarketData[
        contextDict`underlying;spot;contextDict`riskFreeRate;contextDict`dividendYield;contextDict`volatility];
    priceRes:.engine.priceOption[legTrade;mktData;model;fdmConfig];
    greeksRes:.greeks.calculateGreeks[legTrade;mktData;model;fdmConfig];
    `unitPrice`unitDelta`unitGamma`unitTheta`remainingTime`isExpiring!(
        priceRes`unitPrice;first greeksRes`delta;first greeksRes`gamma;first greeksRes`theta;rt;isExp)
 };

/ Mark every leg in legs at currentTimeYears. Returns a table indexed in row order
/ matching legs (joined back via legId in step).
.strategy.calendarRoll.__markLegs:{[legs;currentTimeYears;contextDict;model;fdmConfig]
    if[0=count legs; :([] legId:0#0; unitPrice:0#0f; unitDelta:0#0f; unitGamma:0#0f; unitTheta:0#0f; remainingTime:0#0f; isExpiring:0#0b)];
    pricer:.strategy.calendarRoll.__priceLeg[;currentTimeYears;contextDict;model;fdmConfig];
    rowDicts:pricer each legs;
    markTbl:.strategy.__rowDictsToTable rowDicts;
    update legId:legs`legId from markTbl
 };

/ Create one new leg, pricing it at the entry time. legSpec must contain
/ legId, legRole, optionType, strike, side, units, expiryTimeYears, entryTimeYears.
.strategy.calendarRoll.__createLeg:{[legSpec;contextDict;model;fdmConfig]
    fullSpec:legSpec,`entrySpot`entryPrice`currentMark`currentDelta`currentGamma`currentTheta!(
        contextDict`spot;0f;0f;0f;0f;0f);
    pricing:.strategy.calendarRoll.__priceLeg[fullSpec;legSpec`entryTimeYears;contextDict;model;fdmConfig];
    @[fullSpec;`entryPrice`currentMark`currentDelta`currentGamma`currentTheta;:;
        (pricing`unitPrice;pricing`unitPrice;pricing`unitDelta;pricing`unitGamma;pricing`unitTheta)]
 };

/ Compute net position metrics (value, delta, gamma, theta) from a leg table.
.strategy.calendarRoll.__netPosition:{[legs]
    if[0=count legs; :`value`delta`gamma`theta!(0f;0f;0f;0f)];
    `value`delta`gamma`theta!(
        sum (legs`side)*(legs`units)*legs`currentMark;
        sum (legs`side)*(legs`units)*legs`currentDelta;
        sum (legs`side)*(legs`units)*legs`currentGamma;
        sum (legs`side)*(legs`units)*legs`currentTheta)
 };

.strategy.calendarRoll.init:{[trade;firstStep;model;fdmConfig;stratCfg]
    .strategy.calendarRoll.__validateConfig stratCfg;
    spot:firstStep`spot;
    vol:firstStep`volatility;
    initialTime:0f;
    units:trade`notional;
    contextDict:`spot`volatility`riskFreeRate`dividendYield`underlying`rollThresholdYears!(
        spot;vol;firstStep`riskFreeRate;firstStep`dividendYield;trade`underlying;stratCfg`rollThresholdYears);
    longCal:(stratCfg`spreadType)=`longCalendar;
    frontSide:$[longCal;-1f;1f];
    backSide:$[longCal;1f;-1f];
    optionTypeVal:stratCfg`optionType;
    strikeVal:trade`strike;
    frontSpec:`legId`legRole`optionType`strike`side`units`expiryTimeYears`entryTimeYears!(
        0;`front;optionTypeVal;strikeVal;frontSide;units;initialTime+stratCfg`frontTenorYears;initialTime);
    backSpec:`legId`legRole`optionType`strike`side`units`expiryTimeYears`entryTimeYears!(
        1;`back;optionTypeVal;strikeVal;backSide;units;initialTime+stratCfg`backTenorYears;initialTime);
    frontLeg:.strategy.calendarRoll.__createLeg[frontSpec;contextDict;model;fdmConfig];
    backLeg:.strategy.calendarRoll.__createLeg[backSpec;contextDict;model;fdmConfig];
    legs:.strategy.__rowDictsToTable (frontLeg;backLeg);
    netPos:.strategy.calendarRoll.__netPosition legs;
    positionValue:netPos`value;
    netDelta:netPos`delta;
    netGamma:netPos`gamma;
    netTheta:netPos`theta;
    legEntryTxnCost:sum (legs`units)*(abs legs`entryPrice)*stratCfg`txnCostRate;
    hedgeDeltaOn:stratCfg`hedgeDelta;
    hedgeInit:$[hedgeDeltaOn;
        .strategy.__hedgeInit `spot`positionDelta`txnCostRate!(spot;netDelta;stratCfg`txnCostRate);
        `hedgePosition`hedgeTrade`txnCost`cashAdj!(0f;0f;0f;0f)];
    initialTxnCost:legEntryTxnCost+hedgeInit`txnCost;
    initialStepPnl:neg initialTxnCost;
    initialCash:((neg positionValue)-legEntryTxnCost)+hedgeInit`cashAdj;
    rowEmit:.strategy.calendarRoll.__rowEmitCols!(
        firstStep`stepIndex;firstStep`stepDate;spot;initialTime;2;
        stratCfg`frontTenorYears;stratCfg`backTenorYears;
        positionValue;netDelta;hedgeInit`hedgePosition;hedgeInit`hedgeTrade;initialTxnCost;
        0f;0f;0f;0f;0f;0;
        initialStepPnl;initialStepPnl;0f;
        `OK;"");
    `legs`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue`prevPositionGamma`prevPositionTheta`currentTimeYears`totalRolls`nextLegId`hedgeTrade`txnCost`financingPnl`hedgePnl`positionPnl`rollPnl`thetaPnl`stepPnl`rollEvents`cumulativePnl`rowEmit!(
        legs;initialCash;hedgeInit`hedgePosition;netDelta;1;spot;positionValue;netGamma;netTheta;initialTime;0;2;
        hedgeInit`hedgeTrade;initialTxnCost;0f;0f;0f;0f;0f;initialStepPnl;0;initialStepPnl;rowEmit)
 };

.strategy.calendarRoll.step:{[state;marketStep;trade;model;fdmConfig;stratCfg]
    spot:marketStep`spot;
    vol:marketStep`volatility;
    stepYears:stratCfg`stepYears;
    rollThr:stratCfg`rollThresholdYears;
    newTimeYears:(state`currentTimeYears)+stepYears;
    contextDict:`spot`volatility`riskFreeRate`dividendYield`underlying`rollThresholdYears!(
        spot;vol;marketStep`riskFreeRate;marketStep`dividendYield;trade`underlying;rollThr);

    oldLegs:state`legs;
    markedRows:.strategy.calendarRoll.__markLegs[oldLegs;newTimeYears;contextDict;model;fdmConfig];
    legsWithNew:oldLegs lj `legId xkey markedRows;
    legsWithNew:update markValueNew:side*units*unitPrice, markValuePrev:side*units*currentMark from legsWithNew;

    expiringLegs:legsWithNew where legsWithNew`isExpiring;
    survivingLegs:legsWithNew where not legsWithNew`isExpiring;

    positionPnl:$[0<count survivingLegs; sum (survivingLegs`markValueNew)-survivingLegs`markValuePrev; 0f];
    rollPnl:$[0<count expiringLegs; sum (expiringLegs`markValueNew)-expiringLegs`markValuePrev; 0f];
    rollEventsVal:count expiringLegs;

    txnCostRate:stratCfg`txnCostRate;
    settlements:$[0<count expiringLegs; sum (expiringLegs`side)*(expiringLegs`units)*expiringLegs`unitPrice; 0f];
    rollTxnCloseCost:$[0<count expiringLegs; sum (expiringLegs`units)*(abs expiringLegs`unitPrice)*txnCostRate; 0f];

    units:trade`notional;
    longCal:(stratCfg`spreadType)=`longCalendar;
    frontSide:$[longCal;-1f;1f];
    backSide:$[longCal;1f;-1f];
    optionTypeVal:stratCfg`optionType;
    strikeVal:trade`strike;
    expiringRoles:$[0<count expiringLegs; expiringLegs`legRole; 0#`];
    nextLegId:state`nextLegId;
    newLegSpecs:();
    if[`front in expiringRoles;
        newLegSpecs,:enlist `legId`legRole`optionType`strike`side`units`expiryTimeYears`entryTimeYears!(
            nextLegId;`front;optionTypeVal;strikeVal;frontSide;units;newTimeYears+stratCfg`frontTenorYears;newTimeYears);
        nextLegId+:1];
    if[(`back in expiringRoles)&stratCfg`rollBackLeg;
        newLegSpecs,:enlist `legId`legRole`optionType`strike`side`units`expiryTimeYears`entryTimeYears!(
            nextLegId;`back;optionTypeVal;strikeVal;backSide;units;newTimeYears+stratCfg`backTenorYears;newTimeYears);
        nextLegId+:1];

    newEntryRows:$[0<count newLegSpecs;
        .strategy.calendarRoll.__createLeg[;contextDict;model;fdmConfig] each newLegSpecs;
        ()];
    newEntryTbl:$[0<count newEntryRows; .strategy.__rowDictsToTable newEntryRows; .strategy.calendarRoll.__emptyLegTable[]];
    newEntryCost:$[0<count newEntryTbl; sum (newEntryTbl`side)*(newEntryTbl`units)*newEntryTbl`entryPrice; 0f];
    rollTxnOpenCost:$[0<count newEntryTbl; sum (newEntryTbl`units)*(abs newEntryTbl`entryPrice)*txnCostRate; 0f];
    rollTxnCost:rollTxnCloseCost+rollTxnOpenCost;
    rollCashFlow:(settlements-newEntryCost)-rollTxnCost;

    totalRolls:(state`totalRolls)+rollEventsVal;

    survivingUpdated:$[0<count survivingLegs;
        select legId,legRole,optionType,strike,side,units,expiryTimeYears,entryTimeYears,entrySpot,entryPrice,
            currentMark:unitPrice,currentDelta:unitDelta,currentGamma:unitGamma,currentTheta:unitTheta
            from survivingLegs;
        .strategy.calendarRoll.__emptyLegTable[]];
    newLegs:$[0<count newEntryTbl; survivingUpdated,newEntryTbl; survivingUpdated];

    netPosNew:.strategy.calendarRoll.__netPosition newLegs;
    newPositionValue:netPosNew`value;
    netDelta:netPosNew`delta;
    netGamma:netPosNew`gamma;
    netTheta:netPosNew`theta;

    cashPrev:state`cash;
    prevHedgePos:state`hedgePosition;
    hedgeDeltaOn:stratCfg`hedgeDelta;
    financingRate:stratCfg`financingRate;

    hedgeUpdate:$[hedgeDeltaOn;
        .strategy.__hedgeStep[
            `cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue!(
                cashPrev;prevHedgePos;state`hedgedDelta;state`numRebalances;state`prevSpot;state`prevPositionValue);
            `spot`positionValue`positionDelta`stepIndex`stepYears`txnCostRate`financingRate`rebalanceMode`rebalanceInterval`deltaBand!(
                spot;newPositionValue;netDelta;marketStep`stepIndex;stepYears;txnCostRate;financingRate;stratCfg`rebalanceMode;stratCfg`rebalanceInterval;stratCfg`deltaBand)];
        `cash`hedgePosition`hedgedDelta`numRebalances`hedgeTrade`txnCost`financingPnl`hedgePnl!(
            cashPrev+(financingRate*cashPrev)*stepYears;
            0f;0f;state`numRebalances;0f;0f;(financingRate*cashPrev)*stepYears;0f)];

    hedgeTradeVal:hedgeUpdate`hedgeTrade;
    hedgeTxnCost:hedgeUpdate`txnCost;
    financingPnl:hedgeUpdate`financingPnl;
    hedgePnl:hedgeUpdate`hedgePnl;
    newHedgePos:hedgeUpdate`hedgePosition;
    newHedgedDelta:hedgeUpdate`hedgedDelta;
    numRebalances:hedgeUpdate`numRebalances;
    cashAfterHedge:hedgeUpdate`cash;
    newCash:cashAfterHedge+rollCashFlow;

    totalTxnCost:hedgeTxnCost+rollTxnCost;
    stepPnl:(positionPnl+rollPnl+hedgePnl+financingPnl)-totalTxnCost;
    spotMove:spot-state`prevSpot;
    theoreticalGammaPnl:(0.5*state`prevPositionGamma)*spotMove*spotMove;
    thetaPnl:(state`prevPositionTheta)*stepYears;
    cumulativePnl:(state`cumulativePnl)+stepPnl;

    activeLegCount:count newLegs;
    sortedByRemaining:$[0<activeLegCount; `expiryTimeYears xasc newLegs; ()];
    frontRemainingTime:$[0<activeLegCount; (first sortedByRemaining)[`expiryTimeYears]-newTimeYears; 0Nf];
    backRemainingTime:$[1<activeLegCount; (sortedByRemaining 1)[`expiryTimeYears]-newTimeYears; 0Nf];

    rowEmit:.strategy.calendarRoll.__rowEmitCols!(
        marketStep`stepIndex;marketStep`stepDate;spot;newTimeYears;activeLegCount;frontRemainingTime;backRemainingTime;
        newPositionValue;netDelta;newHedgePos;hedgeTradeVal;totalTxnCost;
        positionPnl;rollPnl;hedgePnl;financingPnl;thetaPnl;rollEventsVal;
        stepPnl;cumulativePnl;theoreticalGammaPnl;
        `OK;"");

    @[state;
        `legs`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue`prevPositionGamma`prevPositionTheta`currentTimeYears`totalRolls`nextLegId`hedgeTrade`txnCost`financingPnl`hedgePnl`positionPnl`rollPnl`thetaPnl`stepPnl`rollEvents`cumulativePnl`rowEmit;:;
        (newLegs;newCash;newHedgePos;newHedgedDelta;numRebalances;spot;newPositionValue;netGamma;netTheta;newTimeYears;totalRolls;nextLegId;hedgeTradeVal;totalTxnCost;financingPnl;hedgePnl;positionPnl;rollPnl;thetaPnl;stepPnl;rollEventsVal;cumulativePnl;rowEmit)]
 };

.strategy.calendarRoll.summary:{[resultTable;stratCfg]
    emptyResult:`strategyName`spreadType`steps`totalRolls`totalPnl`positionPnlTotal`rollPnlTotal`hedgePnlTotal`financingTotal`txnCostTotal`numRebalances`theoreticalGammaPnlTotal`thetaPnlTotal`gammaReconResidual`maxDrawdown`meanStepPnl`stepPnlVol`status`errorMessage!(
        `calendarRoll;stratCfg`spreadType;0;0;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;`ERROR;"empty result table");
    if[(0=count resultTable)|not 98h=type resultTable; :emptyResult];
    okMask:(resultTable`status)=`OK;
    okRows:resultTable where okMask;
    stepCount:count resultTable;
    if[0=count okRows; :@[emptyResult;`steps;:;stepCount]];

    totalsDict:first 0!select
        totalPnl:sum stepPnl,
        positionPnlTotal:sum positionPnl,
        rollPnlTotal:sum rollPnl,
        hedgePnlTotal:sum hedgePnl,
        financingTotal:sum financingPnl,
        txnCostTotal:sum txnCost,
        theoreticalGammaPnlTotal:sum theoreticalGammaPnl,
        thetaPnlTotal:sum thetaPnl,
        meanStepPnl:avg stepPnl,
        stepPnlVol:dev stepPnl
        from okRows;

    nonRollRows:okRows where 0=okRows`rollEvents;
    reconDict:first 0!select
        nonRollStepPnl:sum stepPnl,
        nonRollTxnCost:sum txnCost,
        nonRollFinancing:sum financingPnl,
        nonRollGamma:sum theoreticalGammaPnl,
        nonRollTheta:sum thetaPnl
        from nonRollRows;
    pnlExclCosts:(reconDict[`nonRollStepPnl]+reconDict`nonRollTxnCost)-reconDict`nonRollFinancing;
    gammaReconResidual:pnlExclCosts-(reconDict[`nonRollGamma]+reconDict`nonRollTheta);

    cumPnlSeries:sums okRows`stepPnl;
    drawdownSeries:(maxs cumPnlSeries)-cumPnlSeries;
    maxDrawdownVal:max drawdownSeries;
    numRebalancesVal:sum 0<>okRows`hedgeTrade;
    totalRollsVal:last okRows[`rollEvents]+'sums okRows`rollEvents;
    totalRollsFromTable:sum okRows`rollEvents;

    totalsDict,`strategyName`spreadType`steps`totalRolls`numRebalances`gammaReconResidual`maxDrawdown`status`errorMessage!(
        `calendarRoll;stratCfg`spreadType;stepCount;totalRollsFromTable;numRebalancesVal;gammaReconResidual;maxDrawdownVal;`OK;"")
 };

.strategy.register[
    `calendarRoll;
    .strategy.calendarRoll.init;
    .strategy.calendarRoll.step;
    .strategy.calendarRoll.summary;
    .strategy.calendarRoll.defaultConfig];

/ ==================================================================
/ 7. Concrete strategy: risk reversal (skew speculation)
/ ==================================================================
/ Two-leg, same expiry, different strikes: long one wing + short the other.
/ Per-leg vol = atmVol + skewSlope * moneynessOffset (linear skew). Trade only when
/ the market skewSlope deviates from a config fairSkew by more than skewMargin.
/ Direction `auto picks the side based on skew sign vs fair; `longCallWing /
/ `longPutWing force the side. Hedged delta via __hedgeStep. Accounting via the
/ portfolio-value identity (PV = cash + signed leg marks + hedgePosition*spot).

.strategy.riskReversal.__rowEmitCols:`stepIndex`stepDate`spot`callStrike`putStrike`callVol`putVol`callPrice`putPrice`positionValue`netDelta`hedgePosition`hedgeTrade`txnCost`positionPnl`hedgePnl`financingPnl`thetaPnl`stepPnl`cumulativePnl`theoreticalGammaPnl`gateOpen`status`message;

.strategy.riskReversal.defaultConfig:{[]
    `riskReversalDirection`wingOffsetPct`skewSlope`fairSkew`skewMargin`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears`hedgeDelta!(
        `auto;0.05;-0.5;-0.2;0.05;`interval;1;0.05;0f;0f;1f%252f;1b)
 };

.strategy.riskReversal.__validateConfig:{[stratCfg]
    requiredKeys:`riskReversalDirection`wingOffsetPct`skewSlope`fairSkew`skewMargin`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears`hedgeDelta;
    missing:requiredKeys where not requiredKeys in key stratCfg;
    if[0<count missing; '"riskReversal config missing keys: ",", " sv string missing];
    if[not stratCfg[`riskReversalDirection] in `auto`longCallWing`longPutWing;
        '"riskReversal direction must be auto, longCallWing, or longPutWing"];
    if[(stratCfg`wingOffsetPct)<=0f; '"riskReversal wingOffsetPct must be positive"];
    if[(stratCfg`skewMargin)<0f; '"riskReversal skewMargin must be non-negative"];
    if[(stratCfg`stepYears)<=0f; '"riskReversal stepYears must be positive"];
 };

/ Decide gateOpen + the actual direction taken given config.
.strategy.riskReversal.__resolveDirection:{[stratCfg]
    deviation:(stratCfg`skewSlope)-stratCfg`fairSkew;
    gateOpen:(abs deviation)>stratCfg`skewMargin;
    direction:stratCfg`riskReversalDirection;
    if[direction=`auto;
        direction:$[deviation<0f;`longCallWing;`longPutWing]];
    `gateOpen`direction!(gateOpen;direction)
 };

/ Build the legs (call wing + put wing) for a given direction. Returns a list of two
/ dicts: legId, legRole, optionType, strike, side, units.
.strategy.riskReversal.__buildLegSpecs:{[direction;wingOffsetPct;spot;notional]
    callStrike:spot*1f+wingOffsetPct;
    putStrike:spot*1f-wingOffsetPct;
    longCall:direction=`longCallWing;
    (
        `legId`legRole`optionType`strike`side`units!(0;`callWing;`call;callStrike;$[longCall;1f;-1f];notional);
        `legId`legRole`optionType`strike`side`units!(1;`putWing;`put;putStrike;$[longCall;-1f;1f];notional)
        )
 };

/ Price one leg (call or put) at the per-strike skew-adjusted vol.
.strategy.riskReversal.__priceLeg:{[legSpec;contextDict;model;fdmConfig]
    spot:contextDict`spot;
    atmVol:contextDict`volatility;
    skewSlope:contextDict`skewSlope;
    expiry:contextDict`expiry;
    rfr:contextDict`riskFreeRate;
    divY:contextDict`dividendYield;
    moneynessOffset:((legSpec`strike)-spot)%spot;
    legVol:atmVol+skewSlope*moneynessOffset;
    legVol:0.0001|legVol;
    legTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
        `rrleg;contextDict`underlying;`equityOption;`european;legSpec`optionType;legSpec`strike;expiry;1f);
    mktData:.market.createFlatMarketData[contextDict`underlying;spot;rfr;divY;legVol];
    priceRes:.engine.priceOption[legTrade;mktData;model;fdmConfig];
    greeksRes:.greeks.calculateGreeks[legTrade;mktData;model;fdmConfig];
    `legVol`unitPrice`unitDelta`unitGamma`unitTheta!(legVol;priceRes`unitPrice;first greeksRes`delta;first greeksRes`gamma;first greeksRes`theta)
 };

.strategy.riskReversal.__flatRow:{[marketStep]
    .strategy.riskReversal.__rowEmitCols!(
        marketStep`stepIndex;marketStep`stepDate;marketStep`spot;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;
        0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;0b;`flat;"Gate closed")
 };

.strategy.riskReversal.init:{[trade;firstStep;model;fdmConfig;stratCfg]
    .strategy.riskReversal.__validateConfig stratCfg;
    notional:trade`notional;
    spot:firstStep`spot;
    expiry:trade`expiry;
    direction:.strategy.riskReversal.__resolveDirection stratCfg;
    gateOpen:direction`gateOpen;
    chosenDirection:direction`direction;
    legSpecs:.strategy.riskReversal.__buildLegSpecs[chosenDirection;stratCfg`wingOffsetPct;spot;notional];
    callSpec:legSpecs 0;
    putSpec:legSpecs 1;
    callStrike:callSpec`strike;
    putStrike:putSpec`strike;
    contextDict:`spot`volatility`skewSlope`expiry`riskFreeRate`dividendYield`underlying!(
        spot;firstStep`volatility;stratCfg`skewSlope;expiry;firstStep`riskFreeRate;firstStep`dividendYield;trade`underlying);
    if[not gateOpen;
        flatRow:.strategy.riskReversal.__flatRow firstStep;
        flatRow:@[flatRow;(`callStrike;`putStrike);:;(callStrike;putStrike)];
        :`gateOpen`direction`notional`callStrike`putStrike`callSide`putSide`expiryYears`impliedAtmVolAtEntry`skewSlope`fairSkew`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue`prevPositionGamma`prevPositionTheta`hedgeTrade`txnCost`financingPnl`hedgePnl`positionPnl`thetaPnl`stepPnl`rollEvents`cumulativePnl`rowEmit!(
            0b;chosenDirection;notional;callStrike;putStrike;callSpec`side;putSpec`side;expiry;firstStep`volatility;stratCfg`skewSlope;stratCfg`fairSkew;
            0f;0f;0f;0;spot;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;0;0f;flatRow)];
    callMark:.strategy.riskReversal.__priceLeg[callSpec;contextDict;model;fdmConfig];
    putMark:.strategy.riskReversal.__priceLeg[putSpec;contextDict;model;fdmConfig];
    callSide:callSpec`side;
    putSide:putSpec`side;
    callPrice:callMark`unitPrice;
    putPrice:putMark`unitPrice;
    positionValue:notional*(callSide*callPrice)+putSide*putPrice;
    netDelta:notional*(callSide*callMark`unitDelta)+putSide*putMark`unitDelta;
    netGamma:notional*(callSide*callMark`unitGamma)+putSide*putMark`unitGamma;
    netTheta:notional*(callSide*callMark`unitTheta)+putSide*putMark`unitTheta;
    legEntryTxnCost:notional*((abs callPrice)+abs putPrice)*stratCfg`txnCostRate;
    hedgeInit:$[stratCfg`hedgeDelta;
        .strategy.__hedgeInit `spot`positionDelta`txnCostRate!(spot;netDelta;stratCfg`txnCostRate);
        `hedgePosition`hedgeTrade`txnCost`cashAdj!(0f;0f;0f;0f)];
    initialTxnCost:legEntryTxnCost+hedgeInit`txnCost;
    initialStepPnl:neg initialTxnCost;
    initialCash:((neg positionValue)-legEntryTxnCost)+hedgeInit`cashAdj;
    rowEmit:.strategy.riskReversal.__rowEmitCols!(
        firstStep`stepIndex;firstStep`stepDate;spot;callStrike;putStrike;callMark`legVol;putMark`legVol;callPrice;putPrice;
        positionValue;netDelta;hedgeInit`hedgePosition;hedgeInit`hedgeTrade;initialTxnCost;
        0f;0f;0f;0f;initialStepPnl;initialStepPnl;0f;1b;`OK;"");
    `gateOpen`direction`notional`callStrike`putStrike`callSide`putSide`expiryYears`impliedAtmVolAtEntry`skewSlope`fairSkew`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue`prevPositionGamma`prevPositionTheta`hedgeTrade`txnCost`financingPnl`hedgePnl`positionPnl`thetaPnl`stepPnl`rollEvents`cumulativePnl`rowEmit!(
        1b;chosenDirection;notional;callStrike;putStrike;callSide;putSide;expiry;firstStep`volatility;stratCfg`skewSlope;stratCfg`fairSkew;
        initialCash;hedgeInit`hedgePosition;netDelta;1;spot;positionValue;netGamma;netTheta;
        hedgeInit`hedgeTrade;initialTxnCost;0f;0f;0f;0f;initialStepPnl;0;initialStepPnl;rowEmit)
 };

.strategy.riskReversal.step:{[state;marketStep;trade;model;fdmConfig;stratCfg]
    if[not state`gateOpen;
        flatRow:.strategy.riskReversal.__flatRow marketStep;
        flatRow:@[flatRow;(`callStrike;`putStrike);:;(state`callStrike;state`putStrike)];
        :@[state;`rowEmit;:;flatRow]];
    notional:state`notional;
    callStrike:state`callStrike;
    putStrike:state`putStrike;
    callSide:state`callSide;
    putSide:state`putSide;
    expiry:state`expiryYears;
    spot:marketStep`spot;
    contextDict:`spot`volatility`skewSlope`expiry`riskFreeRate`dividendYield`underlying!(
        spot;marketStep`volatility;stratCfg`skewSlope;expiry;marketStep`riskFreeRate;marketStep`dividendYield;trade`underlying);
    callSpec:`legId`legRole`optionType`strike`side`units!(0;`callWing;`call;callStrike;callSide;notional);
    putSpec:`legId`legRole`optionType`strike`side`units!(1;`putWing;`put;putStrike;putSide;notional);
    callMark:.strategy.riskReversal.__priceLeg[callSpec;contextDict;model;fdmConfig];
    putMark:.strategy.riskReversal.__priceLeg[putSpec;contextDict;model;fdmConfig];
    callPrice:callMark`unitPrice;
    putPrice:putMark`unitPrice;
    newPositionValue:notional*(callSide*callPrice)+putSide*putPrice;
    netDelta:notional*(callSide*callMark`unitDelta)+putSide*putMark`unitDelta;
    netGamma:notional*(callSide*callMark`unitGamma)+putSide*putMark`unitGamma;
    netTheta:notional*(callSide*callMark`unitTheta)+putSide*putMark`unitTheta;
    cashPrev:state`cash;
    financingRate:stratCfg`financingRate;
    hedgeUpdate:$[stratCfg`hedgeDelta;
        .strategy.__hedgeStep[
            `cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue!(
                cashPrev;state`hedgePosition;state`hedgedDelta;state`numRebalances;state`prevSpot;state`prevPositionValue);
            `spot`positionValue`positionDelta`stepIndex`stepYears`txnCostRate`financingRate`rebalanceMode`rebalanceInterval`deltaBand!(
                spot;newPositionValue;netDelta;marketStep`stepIndex;stratCfg`stepYears;stratCfg`txnCostRate;financingRate;stratCfg`rebalanceMode;stratCfg`rebalanceInterval;stratCfg`deltaBand)];
        `cash`hedgePosition`hedgedDelta`numRebalances`hedgeTrade`txnCost`financingPnl`hedgePnl!(
            cashPrev+(financingRate*cashPrev)*stratCfg`stepYears;0f;0f;state`numRebalances;0f;0f;(financingRate*cashPrev)*stratCfg`stepYears;0f)];
    positionPnl:newPositionValue-state`prevPositionValue;
    spotMove:spot-state`prevSpot;
    theoreticalGammaPnl:(0.5*state`prevPositionGamma)*spotMove*spotMove;
    thetaPnl:(state`prevPositionTheta)*stratCfg`stepYears;
    txnCostVal:hedgeUpdate`txnCost;
    stepPnl:(positionPnl+(hedgeUpdate`hedgePnl)+hedgeUpdate`financingPnl)-txnCostVal;
    cumulativePnl:(state`cumulativePnl)+stepPnl;
    rowEmit:.strategy.riskReversal.__rowEmitCols!(
        marketStep`stepIndex;marketStep`stepDate;spot;callStrike;putStrike;callMark`legVol;putMark`legVol;callPrice;putPrice;
        newPositionValue;netDelta;hedgeUpdate`hedgePosition;hedgeUpdate`hedgeTrade;txnCostVal;
        positionPnl;hedgeUpdate`hedgePnl;hedgeUpdate`financingPnl;thetaPnl;stepPnl;cumulativePnl;theoreticalGammaPnl;1b;`OK;"");
    @[state;
        `cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue`prevPositionGamma`prevPositionTheta`hedgeTrade`txnCost`financingPnl`hedgePnl`positionPnl`thetaPnl`stepPnl`cumulativePnl`rowEmit;:;
        (hedgeUpdate`cash;hedgeUpdate`hedgePosition;hedgeUpdate`hedgedDelta;hedgeUpdate`numRebalances;spot;newPositionValue;netGamma;netTheta;
         hedgeUpdate`hedgeTrade;txnCostVal;hedgeUpdate`financingPnl;hedgeUpdate`hedgePnl;positionPnl;thetaPnl;stepPnl;cumulativePnl;rowEmit)]
 };

.strategy.riskReversal.summary:{[resultTable;stratCfg]
    base:`strategyName`steps`gateOpen`direction`skewSlope`fairSkew`callWingVol`putWingVol`totalPnl`positionPnlTotal`hedgePnlTotal`financingTotal`txnCostTotal`numRebalances`theoreticalGammaPnlTotal`thetaPnlTotal`gammaReconResidual`maxDrawdown`meanStepPnl`stepPnlVol`status`errorMessage!(
        `riskReversal;0;0b;`auto;stratCfg`skewSlope;stratCfg`fairSkew;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;`ERROR;"empty result");
    if[(0=count resultTable)|not 98h=type resultTable; :base];
    statusCol:resultTable`status;
    stepCount:count resultTable;
    if[all statusCol=`flat;
        :@[base;(`gateOpen;`steps;`status;`errorMessage;`callWingVol;`putWingVol;`totalPnl;`positionPnlTotal;`hedgePnlTotal;`financingTotal;`txnCostTotal;`theoreticalGammaPnlTotal;`thetaPnlTotal;`gammaReconResidual;`maxDrawdown;`meanStepPnl;`stepPnlVol);:;(0b;stepCount;`flat;"Gate closed";first resultTable`callVol;first resultTable`putVol;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f)]];
    okRows:resultTable where statusCol=`OK;
    if[0=count okRows; :@[base;`steps;:;stepCount]];
    totalsDict:first 0!select
        totalPnl:sum stepPnl,
        positionPnlTotal:sum positionPnl,
        hedgePnlTotal:sum hedgePnl,
        financingTotal:sum financingPnl,
        txnCostTotal:sum txnCost,
        theoreticalGammaPnlTotal:sum theoreticalGammaPnl,
        thetaPnlTotal:sum thetaPnl,
        meanStepPnl:avg stepPnl,
        stepPnlVol:dev stepPnl
        from okRows;
    pnlExclCosts:(totalsDict[`totalPnl]+totalsDict`txnCostTotal)-totalsDict`financingTotal;
    gammaReconResidual:pnlExclCosts-(totalsDict[`theoreticalGammaPnlTotal]+totalsDict`thetaPnlTotal);
    cumPnlSeries:sums okRows`stepPnl;
    maxDrawdownVal:max (maxs cumPnlSeries)-cumPnlSeries;
    numRebalancesVal:sum 0<>okRows`hedgeTrade;
    totalsDict,`strategyName`steps`gateOpen`direction`skewSlope`fairSkew`callWingVol`putWingVol`numRebalances`gammaReconResidual`maxDrawdown`status`errorMessage!(
        `riskReversal;stepCount;1b;`auto;stratCfg`skewSlope;stratCfg`fairSkew;first okRows`callVol;first okRows`putVol;numRebalancesVal;gammaReconResidual;maxDrawdownVal;`OK;"")
 };

.strategy.register[
    `riskReversal;
    .strategy.riskReversal.init;
    .strategy.riskReversal.step;
    .strategy.riskReversal.summary;
    .strategy.riskReversal.defaultConfig];

/ ==================================================================
/ 8. Concrete strategy: model disagreement (two-model entry signal)
/ ==================================================================
/ At init only, price the trade under two model configurations (BS at vol +
/ modelAVolBump vs BS at vol + modelBVolBump). disagreement = priceB - priceA.
/ Gate open when |disagreement| > disagreementThreshold; direction=auto goes long
/ if modelB > modelA, short otherwise. Mark the held position at modelA throughout
/ the path (no per-step MC). Hedge net delta via __hedgeStep. Accounting via
/ portfolio-value identity.
/ Assumption: in this milestone modelA and modelB are both BS with different vol
/ bumps so we can reuse the engine; a true two-model setup (e.g. heston vs BS)
/ would slot in via the same hook with a different per-leg pricer.

.strategy.modelDisagreement.__rowEmitCols:`stepIndex`stepDate`spot`volatility`optionPrice`delta`positionValue`netDelta`hedgePosition`hedgeTrade`txnCost`positionPnl`hedgePnl`financingPnl`thetaPnl`stepPnl`cumulativePnl`theoreticalGammaPnl`gateOpen`tradeSide`status`message;

.strategy.modelDisagreement.defaultConfig:{[]
    `modelAVolBump`modelBVolBump`disagreementThreshold`direction`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears`hedgeDelta!(
        0f;0.05;0.30;`auto;`interval;1;0.05;0f;0f;1f%252f;1b)
 };

.strategy.modelDisagreement.__validateConfig:{[stratCfg]
    requiredKeys:`modelAVolBump`modelBVolBump`disagreementThreshold`direction`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears`hedgeDelta;
    missing:requiredKeys where not requiredKeys in key stratCfg;
    if[0<count missing; '"modelDisagreement config missing keys: ",", " sv string missing];
    if[not stratCfg[`direction] in `auto`long`short;
        '"modelDisagreement direction must be auto, long, or short"];
    if[(stratCfg`disagreementThreshold)<0f;
        '"modelDisagreement disagreementThreshold must be non-negative"];
    if[(stratCfg`stepYears)<=0f; '"modelDisagreement stepYears must be positive"];
 };

.strategy.modelDisagreement.__priceWithVol:{[trade;contextDict;bumpVol;model;fdmConfig]
    volUsed:0.0001|(contextDict`volatility)+bumpVol;
    mktData:.market.createFlatMarketData[contextDict`underlying;contextDict`spot;contextDict`riskFreeRate;contextDict`dividendYield;volUsed];
    priceRes:.engine.priceOption[trade;mktData;model;fdmConfig];
    greeksRes:.greeks.calculateGreeks[trade;mktData;model;fdmConfig];
    `volUsed`unitPrice`unitDelta`unitGamma`unitTheta!(volUsed;priceRes`unitPrice;first greeksRes`delta;first greeksRes`gamma;first greeksRes`theta)
 };

.strategy.modelDisagreement.__flatRow:{[marketStep]
    .strategy.modelDisagreement.__rowEmitCols!(
        marketStep`stepIndex;marketStep`stepDate;marketStep`spot;marketStep`volatility;0Nf;0Nf;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;0b;0;`flat;"Gate closed")
 };

.strategy.modelDisagreement.init:{[trade;firstStep;model;fdmConfig;stratCfg]
    .strategy.modelDisagreement.__validateConfig stratCfg;
    notional:trade`notional;
    spot:firstStep`spot;
    contextDict:`spot`volatility`riskFreeRate`dividendYield`underlying!(
        spot;firstStep`volatility;firstStep`riskFreeRate;firstStep`dividendYield;trade`underlying);
    markA:.strategy.modelDisagreement.__priceWithVol[trade;contextDict;stratCfg`modelAVolBump;model;fdmConfig];
    markB:.strategy.modelDisagreement.__priceWithVol[trade;contextDict;stratCfg`modelBVolBump;model;fdmConfig];
    priceA0:markA`unitPrice;
    priceB0:markB`unitPrice;
    disagreement:priceB0-priceA0;
    gateOpen:(abs disagreement)>stratCfg`disagreementThreshold;
    directionCfg:stratCfg`direction;
    tradeSide:$[directionCfg=`auto; $[disagreement>0f;1f;-1f]; $[directionCfg=`long;1f;-1f]];
    if[not gateOpen;
        flatRow:.strategy.modelDisagreement.__flatRow firstStep;
        :`gateOpen`tradeSide`notional`priceA0`priceB0`disagreement`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue`prevPositionGamma`prevPositionTheta`hedgeTrade`txnCost`financingPnl`hedgePnl`positionPnl`thetaPnl`stepPnl`cumulativePnl`rowEmit!(
            0b;0f;notional;priceA0;priceB0;disagreement;0f;0f;0f;0;spot;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;flatRow)];
    optionUnits:tradeSide*notional;
    positionValue:optionUnits*priceA0;
    positionDelta:optionUnits*markA`unitDelta;
    positionGamma:optionUnits*markA`unitGamma;
    positionTheta:optionUnits*markA`unitTheta;
    legEntryTxnCost:(abs optionUnits)*priceA0*stratCfg`txnCostRate;
    hedgeInit:$[stratCfg`hedgeDelta;
        .strategy.__hedgeInit `spot`positionDelta`txnCostRate!(spot;positionDelta;stratCfg`txnCostRate);
        `hedgePosition`hedgeTrade`txnCost`cashAdj!(0f;0f;0f;0f)];
    initialTxnCost:legEntryTxnCost+hedgeInit`txnCost;
    initialStepPnl:neg initialTxnCost;
    initialCash:((neg positionValue)-legEntryTxnCost)+hedgeInit`cashAdj;
    rowEmit:.strategy.modelDisagreement.__rowEmitCols!(
        firstStep`stepIndex;firstStep`stepDate;spot;firstStep`volatility;priceA0;markA`unitDelta;
        positionValue;positionDelta;hedgeInit`hedgePosition;hedgeInit`hedgeTrade;initialTxnCost;
        0f;0f;0f;0f;initialStepPnl;initialStepPnl;0f;1b;`long$tradeSide;`OK;"");
    `gateOpen`tradeSide`notional`priceA0`priceB0`disagreement`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue`prevPositionGamma`prevPositionTheta`hedgeTrade`txnCost`financingPnl`hedgePnl`positionPnl`thetaPnl`stepPnl`cumulativePnl`rowEmit!(
        1b;tradeSide;notional;priceA0;priceB0;disagreement;initialCash;hedgeInit`hedgePosition;positionDelta;1;spot;positionValue;positionGamma;positionTheta;
        hedgeInit`hedgeTrade;initialTxnCost;0f;0f;0f;0f;initialStepPnl;initialStepPnl;rowEmit)
 };

.strategy.modelDisagreement.step:{[state;marketStep;trade;model;fdmConfig;stratCfg]
    if[not state`gateOpen;
        flatRow:.strategy.modelDisagreement.__flatRow marketStep;
        :@[state;`rowEmit;:;flatRow]];
    notional:state`notional;
    tradeSide:state`tradeSide;
    optionUnits:tradeSide*notional;
    spot:marketStep`spot;
    contextDict:`spot`volatility`riskFreeRate`dividendYield`underlying!(
        spot;marketStep`volatility;marketStep`riskFreeRate;marketStep`dividendYield;trade`underlying);
    markA:.strategy.modelDisagreement.__priceWithVol[trade;contextDict;stratCfg`modelAVolBump;model;fdmConfig];
    optionPrice:markA`unitPrice;
    deltaVal:markA`unitDelta;
    newPositionValue:optionUnits*optionPrice;
    positionDelta:optionUnits*deltaVal;
    positionGamma:optionUnits*markA`unitGamma;
    positionTheta:optionUnits*markA`unitTheta;
    cashPrev:state`cash;
    financingRate:stratCfg`financingRate;
    hedgeUpdate:$[stratCfg`hedgeDelta;
        .strategy.__hedgeStep[
            `cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue!(
                cashPrev;state`hedgePosition;state`hedgedDelta;state`numRebalances;state`prevSpot;state`prevPositionValue);
            `spot`positionValue`positionDelta`stepIndex`stepYears`txnCostRate`financingRate`rebalanceMode`rebalanceInterval`deltaBand!(
                spot;newPositionValue;positionDelta;marketStep`stepIndex;stratCfg`stepYears;stratCfg`txnCostRate;financingRate;stratCfg`rebalanceMode;stratCfg`rebalanceInterval;stratCfg`deltaBand)];
        `cash`hedgePosition`hedgedDelta`numRebalances`hedgeTrade`txnCost`financingPnl`hedgePnl!(
            cashPrev+(financingRate*cashPrev)*stratCfg`stepYears;0f;0f;state`numRebalances;0f;0f;(financingRate*cashPrev)*stratCfg`stepYears;0f)];
    positionPnl:newPositionValue-state`prevPositionValue;
    spotMove:spot-state`prevSpot;
    theoreticalGammaPnl:(0.5*state`prevPositionGamma)*spotMove*spotMove;
    thetaPnl:(state`prevPositionTheta)*stratCfg`stepYears;
    txnCostVal:hedgeUpdate`txnCost;
    stepPnl:(positionPnl+(hedgeUpdate`hedgePnl)+hedgeUpdate`financingPnl)-txnCostVal;
    cumulativePnl:(state`cumulativePnl)+stepPnl;
    rowEmit:.strategy.modelDisagreement.__rowEmitCols!(
        marketStep`stepIndex;marketStep`stepDate;spot;marketStep`volatility;optionPrice;deltaVal;
        newPositionValue;positionDelta;hedgeUpdate`hedgePosition;hedgeUpdate`hedgeTrade;txnCostVal;
        positionPnl;hedgeUpdate`hedgePnl;hedgeUpdate`financingPnl;thetaPnl;stepPnl;cumulativePnl;theoreticalGammaPnl;1b;`long$tradeSide;`OK;"");
    @[state;
        `cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue`prevPositionGamma`prevPositionTheta`hedgeTrade`txnCost`financingPnl`hedgePnl`positionPnl`thetaPnl`stepPnl`cumulativePnl`rowEmit;:;
        (hedgeUpdate`cash;hedgeUpdate`hedgePosition;hedgeUpdate`hedgedDelta;hedgeUpdate`numRebalances;spot;newPositionValue;positionGamma;positionTheta;
         hedgeUpdate`hedgeTrade;txnCostVal;hedgeUpdate`financingPnl;hedgeUpdate`hedgePnl;positionPnl;thetaPnl;stepPnl;cumulativePnl;rowEmit)]
 };

.strategy.modelDisagreement.summary:{[resultTable;stratCfg]
    base:`strategyName`steps`gateOpen`tradeSide`disagreement`totalPnl`positionPnlTotal`hedgePnlTotal`financingTotal`txnCostTotal`numRebalances`theoreticalGammaPnlTotal`thetaPnlTotal`gammaReconResidual`maxDrawdown`status`errorMessage!(
        `modelDisagreement;0;0b;0;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0;0Nf;0Nf;0Nf;0Nf;`ERROR;"empty result");
    if[(0=count resultTable)|not 98h=type resultTable; :base];
    statusCol:resultTable`status;
    stepCount:count resultTable;
    if[all statusCol=`flat;
        :@[base;(`gateOpen;`steps;`status;`errorMessage;`totalPnl;`positionPnlTotal;`hedgePnlTotal;`financingTotal;`txnCostTotal;`theoreticalGammaPnlTotal;`thetaPnlTotal;`gammaReconResidual;`maxDrawdown);:;(0b;stepCount;`flat;"Gate closed";0f;0f;0f;0f;0f;0f;0f;0f;0f)]];
    okRows:resultTable where statusCol=`OK;
    if[0=count okRows; :@[base;`steps;:;stepCount]];
    totalsDict:first 0!select
        totalPnl:sum stepPnl,
        positionPnlTotal:sum positionPnl,
        hedgePnlTotal:sum hedgePnl,
        financingTotal:sum financingPnl,
        txnCostTotal:sum txnCost,
        theoreticalGammaPnlTotal:sum theoreticalGammaPnl,
        thetaPnlTotal:sum thetaPnl
        from okRows;
    pnlExclCosts:(totalsDict[`totalPnl]+totalsDict`txnCostTotal)-totalsDict`financingTotal;
    gammaReconResidual:pnlExclCosts-(totalsDict[`theoreticalGammaPnlTotal]+totalsDict`thetaPnlTotal);
    cumPnlSeries:sums okRows`stepPnl;
    maxDrawdownVal:max (maxs cumPnlSeries)-cumPnlSeries;
    numRebalancesVal:sum 0<>okRows`hedgeTrade;
    totalsDict,`strategyName`steps`gateOpen`tradeSide`disagreement`numRebalances`gammaReconResidual`maxDrawdown`status`errorMessage!(
        `modelDisagreement;stepCount;1b;first okRows`tradeSide;0Nf;numRebalancesVal;gammaReconResidual;maxDrawdownVal;`OK;"")
 };

.strategy.register[
    `modelDisagreement;
    .strategy.modelDisagreement.init;
    .strategy.modelDisagreement.step;
    .strategy.modelDisagreement.summary;
    .strategy.modelDisagreement.defaultConfig];

/ ==================================================================
/ 9. Concrete strategy: delta-vega neutral hedge
/ ==================================================================
/ Hold a primary option (the input trade) and continuously hedge BOTH delta and vega.
/ Hedge instrument: a second vanilla option at config strike (default trade.strike*
/ hedgeStrikeRatio) and config expiry (default trade.expiry). Each step: re-price book
/ and hedge option; vegaHedgeUnits = - bookVega / hedgeOptionVega; netDelta = bookDelta
/ + vegaHedgeUnits*hedgeOptionDelta; delta-hedge netDelta via __hedgeStep. The vega-
/ hedge unit change is a "mini-roll": cash flows -(deltaUnits * hedgeOptionPrice),
/ txnCost on the traded notional. Portfolio-value identity:
/   stepPnl = bookPnl + vegaHedgePnl + hedgePnl + financingPnl - txnCost
/   where txnCost includes both delta-hedge txn cost and vega-rebalance txn cost.

.strategy.deltaVegaHedge.__rowEmitCols:`stepIndex`stepDate`spot`volatility`bookPrice`hedgeOptionPrice`vegaHedgeUnits`positionValue`netDelta`hedgePosition`hedgeTrade`txnCost`bookPnl`vegaHedgePnl`hedgePnl`financingPnl`thetaPnl`stepPnl`cumulativePnl`theoreticalGammaPnl`residualVega`status`message;

.strategy.deltaVegaHedge.defaultConfig:{[]
    `hedgeOptionType`hedgeStrikeRatio`hedgeOptionExpiryYears`vegaRebalanceMode`vegaRebalanceInterval`vegaBand`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears!(
        `call;1.05;0n;`interval;1;0.05;`interval;1;0.05;0f;0f;1f%252f)
 };

.strategy.deltaVegaHedge.__validateConfig:{[stratCfg]
    requiredKeys:`hedgeOptionType`hedgeStrikeRatio`hedgeOptionExpiryYears`vegaRebalanceMode`vegaRebalanceInterval`vegaBand`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears;
    missing:requiredKeys where not requiredKeys in key stratCfg;
    if[0<count missing; '"deltaVegaHedge config missing keys: ",", " sv string missing];
    if[not stratCfg[`hedgeOptionType] in `call`put;
        '"deltaVegaHedge hedgeOptionType must be call or put"];
    if[(stratCfg`hedgeStrikeRatio)<=0f; '"deltaVegaHedge hedgeStrikeRatio must be positive"];
    if[(stratCfg`stepYears)<=0f; '"deltaVegaHedge stepYears must be positive"];
 };

.strategy.deltaVegaHedge.__priceOpt:{[optTrade;contextDict;model;fdmConfig]
    spot:contextDict`spot;
    vol:contextDict`volatility;
    rfr:contextDict`riskFreeRate;
    divY:contextDict`dividendYield;
    mktData:.market.createFlatMarketData[contextDict`underlying;spot;rfr;divY;vol];
    priceRes:.engine.priceOption[optTrade;mktData;model;fdmConfig];
    greeksRes:.greeks.calculateGreeks[optTrade;mktData;model;fdmConfig];
    `unitPrice`unitDelta`unitGamma`unitTheta`unitVega!(priceRes`unitPrice;first greeksRes`delta;first greeksRes`gamma;first greeksRes`theta;first greeksRes`vega)
 };

.strategy.deltaVegaHedge.__buildHedgeOpt:{[bookTrade;stratCfg]
    hedgeStrike:(bookTrade`strike)*stratCfg`hedgeStrikeRatio;
    hedgeExpiry:$[null stratCfg`hedgeOptionExpiryYears; bookTrade`expiry; stratCfg`hedgeOptionExpiryYears];
    @[bookTrade;(`tradeId;`optionType;`strike;`expiry;`notional);:;(`$(string bookTrade`tradeId),"_VH";stratCfg`hedgeOptionType;hedgeStrike;hedgeExpiry;1f)]
 };

.strategy.deltaVegaHedge.init:{[trade;firstStep;model;fdmConfig;stratCfg]
    .strategy.deltaVegaHedge.__validateConfig stratCfg;
    notional:trade`notional;
    spot:firstStep`spot;
    hedgeOpt:.strategy.deltaVegaHedge.__buildHedgeOpt[trade;stratCfg];
    contextDict:`spot`volatility`riskFreeRate`dividendYield`underlying!(
        spot;firstStep`volatility;firstStep`riskFreeRate;firstStep`dividendYield;trade`underlying);
    bookMark:.strategy.deltaVegaHedge.__priceOpt[trade;contextDict;model;fdmConfig];
    hedgeOptMark:.strategy.deltaVegaHedge.__priceOpt[hedgeOpt;contextDict;model;fdmConfig];
    bookPrice0:bookMark`unitPrice;
    hedgeOptPrice0:hedgeOptMark`unitPrice;
    bookValue:notional*bookPrice0;
    bookDelta:notional*bookMark`unitDelta;
    bookVega:notional*bookMark`unitVega;
    bookGamma:notional*bookMark`unitGamma;
    bookTheta:notional*bookMark`unitTheta;
    vegaHedgeUnits:$[(hedgeOptMark`unitVega)=0f; 0f; neg bookVega%hedgeOptMark`unitVega];
    vegaHedgeValue:vegaHedgeUnits*hedgeOptPrice0;
    positionValue:bookValue+vegaHedgeValue;
    netDelta:bookDelta+vegaHedgeUnits*hedgeOptMark`unitDelta;
    netGamma:bookGamma+vegaHedgeUnits*hedgeOptMark`unitGamma;
    netTheta:bookTheta+vegaHedgeUnits*hedgeOptMark`unitTheta;
    residualVega:bookVega+vegaHedgeUnits*hedgeOptMark`unitVega;
    bookEntryTxnCost:(abs bookValue)*stratCfg`txnCostRate;
    vegaEntryTxnCost:(abs vegaHedgeValue)*stratCfg`txnCostRate;
    hedgeInit:.strategy.__hedgeInit `spot`positionDelta`txnCostRate!(spot;netDelta;stratCfg`txnCostRate);
    initialTxnCost:(bookEntryTxnCost+vegaEntryTxnCost)+hedgeInit`txnCost;
    initialStepPnl:neg initialTxnCost;
    initialCash:(((neg positionValue)-bookEntryTxnCost)-vegaEntryTxnCost)+hedgeInit`cashAdj;
    rowEmit:.strategy.deltaVegaHedge.__rowEmitCols!(
        firstStep`stepIndex;firstStep`stepDate;spot;firstStep`volatility;bookPrice0;hedgeOptPrice0;vegaHedgeUnits;
        positionValue;netDelta;hedgeInit`hedgePosition;hedgeInit`hedgeTrade;initialTxnCost;
        0f;0f;0f;0f;0f;initialStepPnl;initialStepPnl;0f;residualVega;`OK;"");
    `bookTrade`hedgeOpt`notional`vegaHedgeUnits`prevBookPrice`prevHedgeOptPrice`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue`prevPositionGamma`prevPositionTheta`hedgeTrade`txnCost`financingPnl`hedgePnl`bookPnl`vegaHedgePnl`thetaPnl`stepPnl`cumulativePnl`rowEmit!(
        trade;hedgeOpt;notional;vegaHedgeUnits;bookPrice0;hedgeOptPrice0;initialCash;hedgeInit`hedgePosition;netDelta;1;spot;positionValue;netGamma;netTheta;
        hedgeInit`hedgeTrade;initialTxnCost;0f;0f;0f;0f;0f;initialStepPnl;initialStepPnl;rowEmit)
 };

.strategy.deltaVegaHedge.step:{[state;marketStep;trade;model;fdmConfig;stratCfg]
    notional:state`notional;
    bookTrade:state`bookTrade;
    hedgeOpt:state`hedgeOpt;
    prevVegaHedgeUnits:state`vegaHedgeUnits;
    prevBookPrice:state`prevBookPrice;
    prevHedgeOptPrice:state`prevHedgeOptPrice;
    spot:marketStep`spot;
    contextDict:`spot`volatility`riskFreeRate`dividendYield`underlying!(
        spot;marketStep`volatility;marketStep`riskFreeRate;marketStep`dividendYield;trade`underlying);
    bookMark:.strategy.deltaVegaHedge.__priceOpt[bookTrade;contextDict;model;fdmConfig];
    hedgeOptMark:.strategy.deltaVegaHedge.__priceOpt[hedgeOpt;contextDict;model;fdmConfig];
    bookPrice:bookMark`unitPrice;
    hedgeOptPrice:hedgeOptMark`unitPrice;
    bookValue:notional*bookPrice;
    bookDelta:notional*bookMark`unitDelta;
    bookVega:notional*bookMark`unitVega;
    bookGamma:notional*bookMark`unitGamma;
    bookTheta:notional*bookMark`unitTheta;
    newVegaHedgeUnits:$[(hedgeOptMark`unitVega)=0f; prevVegaHedgeUnits; neg bookVega%hedgeOptMark`unitVega];
    vegaUnitsChange:newVegaHedgeUnits-prevVegaHedgeUnits;
    bookPnl:notional*(bookPrice-prevBookPrice);
    vegaHedgePnl:prevVegaHedgeUnits*(hedgeOptPrice-prevHedgeOptPrice);
    vegaTradeCash:neg vegaUnitsChange*hedgeOptPrice;
    vegaTradeTxn:(abs vegaUnitsChange*hedgeOptPrice)*stratCfg`txnCostRate;
    newVegaHedgeValue:newVegaHedgeUnits*hedgeOptPrice;
    newPositionValue:bookValue+newVegaHedgeValue;
    netDelta:bookDelta+newVegaHedgeUnits*hedgeOptMark`unitDelta;
    netGamma:bookGamma+newVegaHedgeUnits*hedgeOptMark`unitGamma;
    netTheta:bookTheta+newVegaHedgeUnits*hedgeOptMark`unitTheta;
    residualVega:bookVega+newVegaHedgeUnits*hedgeOptMark`unitVega;
    cashPrev:state`cash;
    financingRate:stratCfg`financingRate;
    hedgeUpdate:.strategy.__hedgeStep[
        `cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue!(
            cashPrev;state`hedgePosition;state`hedgedDelta;state`numRebalances;state`prevSpot;state`prevPositionValue);
        `spot`positionValue`positionDelta`stepIndex`stepYears`txnCostRate`financingRate`rebalanceMode`rebalanceInterval`deltaBand!(
            spot;newPositionValue;netDelta;marketStep`stepIndex;stratCfg`stepYears;stratCfg`txnCostRate;financingRate;stratCfg`rebalanceMode;stratCfg`rebalanceInterval;stratCfg`deltaBand)];
    hedgeTxnCost:hedgeUpdate`txnCost;
    totalTxnCost:hedgeTxnCost+vegaTradeTxn;
    newCash:(hedgeUpdate`cash)+vegaTradeCash-vegaTradeTxn;
    spotMove:spot-state`prevSpot;
    theoreticalGammaPnl:(0.5*state`prevPositionGamma)*spotMove*spotMove;
    thetaPnl:(state`prevPositionTheta)*stratCfg`stepYears;
    stepPnl:(bookPnl+vegaHedgePnl+(hedgeUpdate`hedgePnl)+hedgeUpdate`financingPnl)-totalTxnCost;
    cumulativePnl:(state`cumulativePnl)+stepPnl;
    rowEmit:.strategy.deltaVegaHedge.__rowEmitCols!(
        marketStep`stepIndex;marketStep`stepDate;spot;marketStep`volatility;bookPrice;hedgeOptPrice;newVegaHedgeUnits;
        newPositionValue;netDelta;hedgeUpdate`hedgePosition;hedgeUpdate`hedgeTrade;totalTxnCost;
        bookPnl;vegaHedgePnl;hedgeUpdate`hedgePnl;hedgeUpdate`financingPnl;thetaPnl;stepPnl;cumulativePnl;theoreticalGammaPnl;residualVega;`OK;"");
    @[state;
        `vegaHedgeUnits`prevBookPrice`prevHedgeOptPrice`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue`prevPositionGamma`prevPositionTheta`hedgeTrade`txnCost`financingPnl`hedgePnl`bookPnl`vegaHedgePnl`thetaPnl`stepPnl`cumulativePnl`rowEmit;:;
        (newVegaHedgeUnits;bookPrice;hedgeOptPrice;newCash;hedgeUpdate`hedgePosition;hedgeUpdate`hedgedDelta;hedgeUpdate`numRebalances;spot;newPositionValue;netGamma;netTheta;
         hedgeUpdate`hedgeTrade;totalTxnCost;hedgeUpdate`financingPnl;hedgeUpdate`hedgePnl;bookPnl;vegaHedgePnl;thetaPnl;stepPnl;cumulativePnl;rowEmit)]
 };

.strategy.deltaVegaHedge.summary:{[resultTable;stratCfg]
    base:`strategyName`steps`totalPnl`bookPnlTotal`vegaHedgePnlTotal`hedgePnlTotal`financingTotal`txnCostTotal`numRebalances`numVegaRehedges`maxAbsResidualVega`bookVega0`theoreticalGammaPnlTotal`thetaPnlTotal`gammaReconResidual`maxDrawdown`status`errorMessage!(
        `deltaVegaHedge;0;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0;0;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;`ERROR;"empty result");
    if[(0=count resultTable)|not 98h=type resultTable; :base];
    okRows:resultTable where (resultTable`status)=`OK;
    stepCount:count resultTable;
    if[0=count okRows; :@[base;`steps;:;stepCount]];
    totalsDict:first 0!select
        totalPnl:sum stepPnl,
        bookPnlTotal:sum bookPnl,
        vegaHedgePnlTotal:sum vegaHedgePnl,
        hedgePnlTotal:sum hedgePnl,
        financingTotal:sum financingPnl,
        txnCostTotal:sum txnCost,
        theoreticalGammaPnlTotal:sum theoreticalGammaPnl,
        thetaPnlTotal:sum thetaPnl
        from okRows;
    pnlExclCosts:(totalsDict[`totalPnl]+totalsDict`txnCostTotal)-totalsDict`financingTotal;
    gammaReconResidual:pnlExclCosts-(totalsDict[`theoreticalGammaPnlTotal]+totalsDict`thetaPnlTotal);
    cumPnlSeries:sums okRows`stepPnl;
    maxDrawdownVal:max (maxs cumPnlSeries)-cumPnlSeries;
    numRebalancesVal:sum 0<>okRows`hedgeTrade;
    vegaUnitsSeries:okRows`vegaHedgeUnits;
    numVegaRehedges:sum 0<>(vegaUnitsSeries-prev vegaUnitsSeries);
    maxAbsResidualVega:max abs okRows`residualVega;
    bookVega0:(0Nf,prev okRows`bookPrice)0;
    bookVega0:0Nf;
    totalsDict,`strategyName`steps`numRebalances`numVegaRehedges`maxAbsResidualVega`bookVega0`gammaReconResidual`maxDrawdown`status`errorMessage!(
        `deltaVegaHedge;stepCount;numRebalancesVal;numVegaRehedges;maxAbsResidualVega;bookVega0;gammaReconResidual;maxDrawdownVal;`OK;"")
 };

.strategy.register[
    `deltaVegaHedge;
    .strategy.deltaVegaHedge.init;
    .strategy.deltaVegaHedge.step;
    .strategy.deltaVegaHedge.summary;
    .strategy.deltaVegaHedge.defaultConfig];

/ ==================================================================
/ 10. Ensemble / multi-path portfolio runner (Part C)
/ ==================================================================
/ Run N strategies x M paths with common random numbers (same path ensemble fed to
/ every strategy) so cross-strategy comparisons are meaningful. Aggregate to a per-
/ strategy distribution + a cross-strategy correlation matrix + a book-level rollup.

.strategy.path.ensemble:{[pathCfg;numPaths;baseSeed]
    if[numPaths<=0; '"strategy.path.ensemble numPaths must be positive"];
    seeds:baseSeed+til numPaths;
    {.strategy.path.fromSynthetic[@[x;`seed;:;y]]}[pathCfg;] each seeds
 };

.strategy.portfolio.__runOne:{[strategySpec;trade;model;fdmConfig;pathBundle]
    pathId:pathBundle 0;
    path:pathBundle 1;
    sn:strategySpec 0;
    sc:strategySpec 1;
    bundleResult:.[.strategy.runAndSummarize;(sn;trade;path;model;fdmConfig;sc);{x}];
    if[10h=type bundleResult;
        :`strategyName`pathId`totalPnl`maxDrawdown`numRebalances`status`errorMessage!(
            sn;pathId;0Nf;0Nf;0;`ERROR;bundleResult)];
    sumDict:bundleResult`summary;
    numRebVal:$[`numRebalances in key sumDict; sumDict`numRebalances; 0];
    msgVal:$[`errorMessage in key sumDict; sumDict`errorMessage; ""];
    `strategyName`pathId`totalPnl`maxDrawdown`numRebalances`status`errorMessage!(
        sn;pathId;sumDict`totalPnl;sumDict`maxDrawdown;numRebVal;sumDict`status;msgVal)
 };

.strategy.portfolio.runEnsemble:{[strategySpecs;trade;pathEnsemble;model;fdmConfig]
    if[0=count strategySpecs; '"strategy.portfolio.runEnsemble: empty strategySpecs"];
    if[0=count pathEnsemble; '"strategy.portfolio.runEnsemble: empty pathEnsemble"];
    numPaths:count pathEnsemble;
    allRows:();
    stratIdx:0;
    while[stratIdx<count strategySpecs;
        stratSpec:strategySpecs stratIdx;
        pathIdx:0;
        while[pathIdx<numPaths;
            row:.strategy.portfolio.__runOne[stratSpec;trade;model;fdmConfig;(pathIdx;pathEnsemble pathIdx)];
            allRows,:enlist row;
            pathIdx+:1];
        stratIdx+:1];
    .strategy.__rowDictsToTable allRows
 };

.strategy.portfolio.__percentile:{[vec;pct]
    if[0=count vec; :0Nf];
    sortedVec:asc vec;
    idx:floor (pct%100f)*-1+count sortedVec;
    sortedVec idx
 };

.strategy.portfolio.performanceByStrategy:{[ensembleSummary]
    emptyResult:([] strategyName:0#`; pathCount:0#0; meanPnl:0#0Nf; stdPnl:0#0Nf; minPnl:0#0Nf; maxPnl:0#0Nf; p05Pnl:0#0Nf; p95Pnl:0#0Nf; sharpeLike:0#0Nf; meanMaxDrawdown:0#0Nf; winRate:0#0Nf);
    if[(0=count ensembleSummary)|not 98h=type ensembleSummary; :emptyResult];
    okRows:ensembleSummary where (ensembleSummary`status)=`OK;
    if[0=count okRows; :emptyResult];
    grouped:0! select
        pathCount:count i,
        meanPnl:avg totalPnl,
        stdPnl:dev totalPnl,
        minPnl:min totalPnl,
        maxPnl:max totalPnl,
        meanMaxDrawdown:avg maxDrawdown,
        winRate:`float$(sum totalPnl>0f)%count i,
        totalPnls:totalPnl
        by strategyName from okRows;
    grouped:update p05Pnl:.strategy.portfolio.__percentile[;5] each totalPnls,
                    p95Pnl:.strategy.portfolio.__percentile[;95] each totalPnls,
                    sharpeLike:?[0f=stdPnl;0Nf;meanPnl%stdPnl]
        from grouped;
    select strategyName,pathCount,meanPnl,stdPnl,minPnl,maxPnl,p05Pnl,p95Pnl,sharpeLike,meanMaxDrawdown,winRate from grouped
 };

.strategy.portfolio.strategyCorrelation:{[ensembleSummary]
    emptyMatrix:`names`matrix!(0#`;());
    if[(0=count ensembleSummary)|not 98h=type ensembleSummary; :emptyMatrix];
    okRows:ensembleSummary where (ensembleSummary`status)=`OK;
    if[0=count okRows; :emptyMatrix];
    pivoted:0! select pnlVec:totalPnl by strategyName from `pathId xasc okRows;
    sNames:pivoted`strategyName;
    sPnls:pivoted`pnlVec;
    n:count sNames;
    if[n<=1; :`names`matrix!(sNames;enlist enlist $[n=0;0n;1f])];
    matrix:{[i;sPnlsLocal]
        viLocal:sPnlsLocal i;
        {[viInner;vjInner]
            if[(0=dev viInner)|0=dev vjInner; :$[viInner~vjInner;1f;0n]];
            cor[viInner;vjInner]
            }[viLocal;] each sPnlsLocal
        }[;sPnls] each til n;
    `names`matrix!(sNames;matrix)
 };

.strategy.portfolio.dashboard:{[ensembleSummary]
    perfTbl:.strategy.portfolio.performanceByStrategy ensembleSummary;
    corDict:.strategy.portfolio.strategyCorrelation ensembleSummary;
    bookAgg:`pathCount`meanBookPnl`stdBookPnl!(0;0Nf;0Nf);
    if[(0=count ensembleSummary)|not 98h=type ensembleSummary;
        :`performanceByStrategy`strategyCorrelation`bookAggregate`dashboardStatus`dashboardMessage!(
            perfTbl;corDict;bookAgg;`ERROR;"Empty or invalid ensembleSummary")];
    okRows:ensembleSummary where (ensembleSummary`status)=`OK;
    dashStatus:`OK;
    dashMsg:"OK";
    if[0=count okRows; dashStatus:`ERROR; dashMsg:"No OK ensemble rows"];
    if[0<count okRows;
        bookPerPath:0! select sumPnl:sum totalPnl by pathId from okRows;
        bookPnlVec:bookPerPath`sumPnl;
        bookAgg:`pathCount`meanBookPnl`stdBookPnl!(count bookPnlVec;avg bookPnlVec;dev bookPnlVec)];
    `performanceByStrategy`strategyCorrelation`bookAggregate`dashboardStatus`dashboardMessage!(
        perfTbl;corDict;bookAgg;dashStatus;dashMsg)
 };

/ ==================================================================
/ 11. Concrete strategy: longVol (mirror of shortVariance)
/ ==================================================================
/ Buys an ATM straddle when forecastVol > impliedVolAtEntry + entryMargin. Long
/ position has positive gamma, negative theta. Pays premium up front.

.strategy.longVol.__rowEmitCols:`stepIndex`stepDate`spot`volatility`callPrice`putPrice`positionValue`netDelta`hedgePosition`hedgeTrade`txnCost`positionPnl`hedgePnl`financingPnl`thetaPnl`stepPnl`cumulativePnl`theoreticalGammaPnl`premiumPaid`status`message;

.strategy.longVol.defaultConfig:{[]
    `forecastVol`entryMargin`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears!(
        0.30;0.02;`interval;1;0.05;0f;0f;1f%252f)
 };

.strategy.longVol.__validateConfig:{[stratCfg]
    requiredKeys:`forecastVol`entryMargin`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears;
    missing:requiredKeys where not requiredKeys in key stratCfg;
    if[0<count missing; '"longVol config missing keys: ",", " sv string missing];
    if[(stratCfg`forecastVol)<0f; '"longVol forecastVol must be non-negative"];
    if[(stratCfg`entryMargin)<0f; '"longVol entryMargin must be non-negative"];
    if[not stratCfg[`rebalanceMode] in `interval`band;
        '"longVol rebalanceMode must be interval or band"];
    if[(stratCfg`stepYears)<=0f; '"longVol stepYears must be positive"];
 };

.strategy.longVol.__buildLegs:{[trade]
    callTrade:@[trade;(`tradeId;`optionType);:;(`$(string trade`tradeId),"_C";`call)];
    putTrade:@[trade;(`tradeId;`optionType);:;(`$(string trade`tradeId),"_P";`put)];
    `callTrade`putTrade!(callTrade;putTrade)
 };

.strategy.longVol.__flatRow:{[marketStep]
    .strategy.longVol.__rowEmitCols!(
        marketStep`stepIndex;marketStep`stepDate;marketStep`spot;marketStep`volatility;
        0Nf;0Nf;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;
        `flat;"Entry gate closed")
 };

.strategy.longVol.init:{[trade;firstStep;model;fdmConfig;stratCfg]
    .strategy.longVol.__validateConfig stratCfg;
    notional:trade`notional;
    spot:firstStep`spot;
    vol:firstStep`volatility;
    rfr:firstStep`riskFreeRate;
    divY:firstStep`dividendYield;
    forecastVol:stratCfg`forecastVol;
    entryMargin:stratCfg`entryMargin;
    gateOpen:forecastVol>vol+entryMargin;
    legs:.strategy.longVol.__buildLegs trade;
    callTrade:legs`callTrade;
    putTrade:legs`putTrade;
    if[not gateOpen;
        flatRow:.strategy.longVol.__flatRow firstStep;
        :`gateOpen`notional`callTrade`putTrade`impliedVolAtEntry`forecastVol`premiumPaid`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue`prevPositionGamma`prevPositionTheta`hedgeTrade`txnCost`financingPnl`hedgePnl`positionPnl`thetaPnl`stepPnl`cumulativePnl`rowEmit!(
            0b;notional;callTrade;putTrade;vol;forecastVol;0f;0f;0f;0f;0;spot;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;flatRow)];
    mktData:.market.createFlatMarketData[trade`underlying;spot;rfr;divY;vol];
    callPrice:(.engine.priceOption[callTrade;mktData;model;fdmConfig])`unitPrice;
    putPrice:(.engine.priceOption[putTrade;mktData;model;fdmConfig])`unitPrice;
    callGreeks:.greeks.calculateGreeks[callTrade;mktData;model;fdmConfig];
    putGreeks:.greeks.calculateGreeks[putTrade;mktData;model;fdmConfig];
    callDelta:first callGreeks`delta;
    putDelta:first putGreeks`delta;
    callGamma:first callGreeks`gamma;
    putGamma:first putGreeks`gamma;
    callTheta:first callGreeks`theta;
    putTheta:first putGreeks`theta;
    positionValue:notional*callPrice+putPrice;
    positionDelta:notional*callDelta+putDelta;
    positionGamma:notional*callGamma+putGamma;
    positionTheta:notional*callTheta+putTheta;
    premiumPaid:notional*callPrice+putPrice;
    hedgeInit:.strategy.__hedgeInit `spot`positionDelta`txnCostRate!(spot;positionDelta;stratCfg`txnCostRate);
    initialCash:(neg premiumPaid)+hedgeInit`cashAdj;
    initialStepPnl:neg hedgeInit`txnCost;
    rowEmit:.strategy.longVol.__rowEmitCols!(
        firstStep`stepIndex;firstStep`stepDate;spot;vol;callPrice;putPrice;positionValue;positionDelta;
        hedgeInit`hedgePosition;hedgeInit`hedgeTrade;hedgeInit`txnCost;
        0f;0f;0f;0f;initialStepPnl;initialStepPnl;0f;premiumPaid;`OK;"");
    `gateOpen`notional`callTrade`putTrade`impliedVolAtEntry`forecastVol`premiumPaid`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue`prevPositionGamma`prevPositionTheta`hedgeTrade`txnCost`financingPnl`hedgePnl`positionPnl`thetaPnl`stepPnl`cumulativePnl`rowEmit!(
        1b;notional;callTrade;putTrade;vol;forecastVol;premiumPaid;initialCash;hedgeInit`hedgePosition;positionDelta;1;spot;positionValue;positionGamma;positionTheta;
        hedgeInit`hedgeTrade;hedgeInit`txnCost;0f;0f;0f;0f;initialStepPnl;initialStepPnl;rowEmit)
 };

.strategy.longVol.step:{[state;marketStep;trade;model;fdmConfig;stratCfg]
    if[not state`gateOpen;
        :@[state;`rowEmit;:;.strategy.longVol.__flatRow marketStep]];
    notional:state`notional;
    callTrade:state`callTrade;
    putTrade:state`putTrade;
    spot:marketStep`spot;
    vol:marketStep`volatility;
    mktData:.market.createFlatMarketData[trade`underlying;spot;marketStep`riskFreeRate;marketStep`dividendYield;vol];
    callPrice:(.engine.priceOption[callTrade;mktData;model;fdmConfig])`unitPrice;
    putPrice:(.engine.priceOption[putTrade;mktData;model;fdmConfig])`unitPrice;
    callGreeks:.greeks.calculateGreeks[callTrade;mktData;model;fdmConfig];
    putGreeks:.greeks.calculateGreeks[putTrade;mktData;model;fdmConfig];
    positionValue:notional*callPrice+putPrice;
    positionDelta:notional*(first callGreeks`delta)+first putGreeks`delta;
    positionGamma:notional*(first callGreeks`gamma)+first putGreeks`gamma;
    positionTheta:notional*(first callGreeks`theta)+first putGreeks`theta;
    hedgeState:`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue#state;
    stepInputs:`spot`positionValue`positionDelta`stepIndex`stepYears`txnCostRate`financingRate`rebalanceMode`rebalanceInterval`deltaBand!(
        spot;positionValue;positionDelta;marketStep`stepIndex;stratCfg`stepYears;stratCfg`txnCostRate;stratCfg`financingRate;stratCfg`rebalanceMode;stratCfg`rebalanceInterval;stratCfg`deltaBand);
    newHedge:.strategy.__hedgeStep[hedgeState;stepInputs];
    spotMove:spot-state`prevSpot;
    theoreticalGammaPnl:(0.5*state`prevPositionGamma)*spotMove*spotMove;
    thetaPnl:(state`prevPositionTheta)*stratCfg`stepYears;
    cumulativePnl:(state`cumulativePnl)+newHedge`stepPnl;
    rowEmit:.strategy.longVol.__rowEmitCols!(
        marketStep`stepIndex;marketStep`stepDate;spot;vol;callPrice;putPrice;positionValue;positionDelta;
        newHedge`hedgePosition;newHedge`hedgeTrade;newHedge`txnCost;
        newHedge`positionPnl;newHedge`hedgePnl;newHedge`financingPnl;thetaPnl;
        newHedge`stepPnl;cumulativePnl;theoreticalGammaPnl;
        state`premiumPaid;`OK;"");
    state,newHedge,`prevPositionGamma`prevPositionTheta`cumulativePnl`rowEmit!(
        positionGamma;positionTheta;cumulativePnl;rowEmit)
 };

.strategy.longVol.summary:{[resultTable;stratCfg]
    base:`strategyName`gateOpen`steps`premiumPaid`totalPnl`positionPnlTotal`hedgePnlTotal`txnCostTotal`financingTotal`numRebalances`impliedVolAtEntry`forecastVol`realizedVol`varianceRiskPremium`theoreticalGammaPnlTotal`thetaPnlTotal`gammaReconResidual`maxDrawdown`meanStepPnl`stepPnlVol`status`errorMessage!(
        `longVol;0b;0;0f;0Nf;0Nf;0Nf;0Nf;0Nf;0;0Nf;stratCfg`forecastVol;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;`ERROR;"empty");
    if[(0=count resultTable)|not 98h=type resultTable; :base];
    statusCol:resultTable`status;
    stepCount:count resultTable;
    if[all statusCol=`flat;
        :@[base;(`gateOpen;`steps;`status;`errorMessage;`totalPnl;`positionPnlTotal;`hedgePnlTotal;`txnCostTotal;`financingTotal;`theoreticalGammaPnlTotal;`thetaPnlTotal;`gammaReconResidual;`maxDrawdown;`meanStepPnl;`stepPnlVol;`impliedVolAtEntry;`realizedVol;`varianceRiskPremium);:;(0b;stepCount;`flat;"Gate closed";0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;first resultTable`volatility;0Nf;0Nf)]];
    okRows:resultTable where statusCol=`OK;
    if[0=count okRows; :@[base;`steps;:;stepCount]];
    totalsDict:first 0!select
        totalPnl:sum stepPnl, positionPnlTotal:sum positionPnl, hedgePnlTotal:sum hedgePnl,
        txnCostTotal:sum txnCost, financingTotal:sum financingPnl,
        theoreticalGammaPnlTotal:sum theoreticalGammaPnl, thetaPnlTotal:sum thetaPnl,
        meanStepPnl:avg stepPnl, stepPnlVol:dev stepPnl
        from okRows;
    spots:okRows`spot;
    logRet:$[1<count spots; 1_(log spots)-prev log spots; ()];
    nonNullRet:logRet where not null logRet;
    realizedVol:$[(0<count nonNullRet)&stratCfg[`stepYears]>0f; (dev nonNullRet)%sqrt stratCfg`stepYears; 0Nf];
    impliedVolAtEntry:first okRows`volatility;
    varianceRiskPremium:realizedVol-impliedVolAtEntry;
    pnlExclCosts:(totalsDict[`totalPnl]+totalsDict`txnCostTotal)-totalsDict`financingTotal;
    gammaReconResidual:pnlExclCosts-(totalsDict[`theoreticalGammaPnlTotal]+totalsDict`thetaPnlTotal);
    cumPnlSeries:sums okRows`stepPnl;
    maxDrawdownVal:max (maxs cumPnlSeries)-cumPnlSeries;
    numRebalancesVal:sum 0<>okRows`hedgeTrade;
    premiumPaidVal:first okRows`premiumPaid;
    totalsDict,`strategyName`gateOpen`steps`premiumPaid`numRebalances`impliedVolAtEntry`forecastVol`realizedVol`varianceRiskPremium`gammaReconResidual`maxDrawdown`status`errorMessage!(
        `longVol;1b;stepCount;premiumPaidVal;numRebalancesVal;impliedVolAtEntry;stratCfg`forecastVol;realizedVol;varianceRiskPremium;gammaReconResidual;maxDrawdownVal;`OK;"")
 };

.strategy.register[`longVol;.strategy.longVol.init;.strategy.longVol.step;.strategy.longVol.summary;.strategy.longVol.defaultConfig];

/ ==================================================================
/ 12. Concrete strategy: collarTailHedge (collar OR tailHedge mode)
/ ==================================================================
/ collar: long underlying (notional units) + long OTM put + short OTM call. Held
/ outright (no external delta hedge); net delta tracked for reporting.
/ tailHedge: long OTM puts only, sized to premiumBudgetPct * notional * spot value.
/ Cap (collar): fixed call strike at spot*(1+callStrikePct); NOT solved for zero-cost.

.strategy.collarTailHedge.__rowEmitCols:`stepIndex`stepDate`spot`mode`underlyingValue`putValue`callValue`positionValue`netDelta`txnCost`positionPnl`financingPnl`thetaPnl`stepPnl`cumulativePnl`theoreticalGammaPnl`putUnits`status`message;

.strategy.collarTailHedge.defaultConfig:{[]
    `mode`putStrikePct`callStrikePct`premiumBudgetPct`txnCostRate`financingRate`stepYears!(
        `collar;0.05;0.05;0.005;0f;0f;1f%252f)
 };

.strategy.collarTailHedge.__validateConfig:{[stratCfg]
    requiredKeys:`mode`putStrikePct`callStrikePct`premiumBudgetPct`txnCostRate`financingRate`stepYears;
    missing:requiredKeys where not requiredKeys in key stratCfg;
    if[0<count missing; '"collarTailHedge config missing keys: ",", " sv string missing];
    if[not stratCfg[`mode] in `collar`tailHedge;
        '"collarTailHedge mode must be collar or tailHedge"];
    if[(stratCfg`putStrikePct)<=0f; '"collarTailHedge putStrikePct must be positive"];
    if[(stratCfg[`mode]=`collar)&(stratCfg`callStrikePct)<=0f;
        '"collarTailHedge callStrikePct must be positive in collar mode"];
    if[(stratCfg[`mode]=`tailHedge)&(stratCfg`premiumBudgetPct)<=0f;
        '"collarTailHedge premiumBudgetPct must be positive in tailHedge mode"];
    if[(stratCfg`stepYears)<=0f; '"collarTailHedge stepYears must be positive"];
 };

.strategy.collarTailHedge.__priceLeg:{[legTrade;contextDict;model;fdmConfig]
    mktData:.market.createFlatMarketData[contextDict`underlying;contextDict`spot;contextDict`riskFreeRate;contextDict`dividendYield;contextDict`volatility];
    priceRes:.engine.priceOption[legTrade;mktData;model;fdmConfig];
    greeksRes:.greeks.calculateGreeks[legTrade;mktData;model;fdmConfig];
    `unitPrice`unitDelta`unitGamma`unitTheta!(priceRes`unitPrice;first greeksRes`delta;first greeksRes`gamma;first greeksRes`theta)
 };

.strategy.collarTailHedge.init:{[trade;firstStep;model;fdmConfig;stratCfg]
    .strategy.collarTailHedge.__validateConfig stratCfg;
    notional:trade`notional;
    spot:firstStep`spot;
    mode:stratCfg`mode;
    isCollar:mode=`collar;
    putStrike:spot*1f-stratCfg`putStrikePct;
    callStrike:spot*1f+stratCfg`callStrikePct;
    contextDict:`spot`volatility`riskFreeRate`dividendYield`underlying!(
        spot;firstStep`volatility;firstStep`riskFreeRate;firstStep`dividendYield;trade`underlying);
    putTrade:@[trade;(`tradeId;`optionType;`strike);:;(`$(string trade`tradeId),"_P";`put;putStrike)];
    callTrade:@[trade;(`tradeId;`optionType;`strike);:;(`$(string trade`tradeId),"_C";`call;callStrike)];
    putMark:.strategy.collarTailHedge.__priceLeg[putTrade;contextDict;model;fdmConfig];
    callMark:$[isCollar; .strategy.collarTailHedge.__priceLeg[callTrade;contextDict;model;fdmConfig]; `unitPrice`unitDelta`unitGamma`unitTheta!(0f;0f;0f;0f)];
    putUnits:$[isCollar; notional; (stratCfg[`premiumBudgetPct]*notional*spot)%putMark`unitPrice];
    callUnits:$[isCollar; notional; 0f];
    underlyingUnits:$[isCollar; notional; 0f];
    putValue:putUnits*putMark`unitPrice;
    callValue:callUnits*callMark`unitPrice;
    underlyingValue:underlyingUnits*spot;
    positionValue:(underlyingValue+putValue)-callValue;
    netDelta:underlyingUnits+(putUnits*putMark`unitDelta)-callUnits*callMark`unitDelta;
    netGamma:(putUnits*putMark`unitGamma)-callUnits*callMark`unitGamma;
    netTheta:(putUnits*putMark`unitTheta)-callUnits*callMark`unitTheta;
    legTxnCost:(((abs underlyingValue)+abs putValue)+abs callValue)*stratCfg`txnCostRate;
    initialStepPnl:neg legTxnCost;
    initialCash:(neg positionValue)-legTxnCost;
    rowEmit:.strategy.collarTailHedge.__rowEmitCols!(
        firstStep`stepIndex;firstStep`stepDate;spot;mode;underlyingValue;putValue;callValue;positionValue;netDelta;legTxnCost;
        0f;0f;0f;initialStepPnl;initialStepPnl;0f;putUnits;`OK;"");
    `mode`underlyingUnits`putUnits`callUnits`putTrade`callTrade`cash`hedgePosition`prevSpot`prevPositionValue`prevPositionGamma`prevPositionTheta`txnCost`financingPnl`positionPnl`thetaPnl`stepPnl`cumulativePnl`rowEmit!(
        mode;underlyingUnits;putUnits;callUnits;putTrade;callTrade;initialCash;0f;spot;positionValue;netGamma;netTheta;legTxnCost;0f;0f;0f;initialStepPnl;initialStepPnl;rowEmit)
 };

.strategy.collarTailHedge.step:{[state;marketStep;trade;model;fdmConfig;stratCfg]
    mode:state`mode;
    underlyingUnits:state`underlyingUnits;
    putUnits:state`putUnits;
    callUnits:state`callUnits;
    putTrade:state`putTrade;
    callTrade:state`callTrade;
    spot:marketStep`spot;
    contextDict:`spot`volatility`riskFreeRate`dividendYield`underlying!(
        spot;marketStep`volatility;marketStep`riskFreeRate;marketStep`dividendYield;trade`underlying);
    putMark:.strategy.collarTailHedge.__priceLeg[putTrade;contextDict;model;fdmConfig];
    callMark:$[mode=`collar; .strategy.collarTailHedge.__priceLeg[callTrade;contextDict;model;fdmConfig]; `unitPrice`unitDelta`unitGamma`unitTheta!(0f;0f;0f;0f)];
    putValue:putUnits*putMark`unitPrice;
    callValue:callUnits*callMark`unitPrice;
    underlyingValue:underlyingUnits*spot;
    newPositionValue:(underlyingValue+putValue)-callValue;
    netDelta:underlyingUnits+(putUnits*putMark`unitDelta)-callUnits*callMark`unitDelta;
    netGamma:(putUnits*putMark`unitGamma)-callUnits*callMark`unitGamma;
    netTheta:(putUnits*putMark`unitTheta)-callUnits*callMark`unitTheta;
    cashPrev:state`cash;
    financingPnl:(stratCfg`financingRate)*cashPrev*stratCfg`stepYears;
    positionPnl:newPositionValue-state`prevPositionValue;
    spotMove:spot-state`prevSpot;
    theoreticalGammaPnl:(0.5*state`prevPositionGamma)*spotMove*spotMove;
    thetaPnl:(state`prevPositionTheta)*stratCfg`stepYears;
    stepPnl:positionPnl+financingPnl;
    newCash:cashPrev+financingPnl;
    cumulativePnl:(state`cumulativePnl)+stepPnl;
    rowEmit:.strategy.collarTailHedge.__rowEmitCols!(
        marketStep`stepIndex;marketStep`stepDate;spot;mode;underlyingValue;putValue;callValue;newPositionValue;netDelta;0f;
        positionPnl;financingPnl;thetaPnl;stepPnl;cumulativePnl;theoreticalGammaPnl;putUnits;`OK;"");
    @[state;`cash`hedgePosition`prevSpot`prevPositionValue`prevPositionGamma`prevPositionTheta`txnCost`financingPnl`positionPnl`thetaPnl`stepPnl`cumulativePnl`rowEmit;:;
        (newCash;0f;spot;newPositionValue;netGamma;netTheta;0f;financingPnl;positionPnl;thetaPnl;stepPnl;cumulativePnl;rowEmit)]
 };

.strategy.collarTailHedge.summary:{[resultTable;stratCfg]
    base:`strategyName`mode`steps`netPremium`protectionFloor`cap`budgetUsed`totalPnl`positionPnlTotal`financingTotal`txnCostTotal`theoreticalGammaPnlTotal`thetaPnlTotal`maxDrawdown`status`errorMessage!(
        `collarTailHedge;stratCfg`mode;0;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;`ERROR;"empty");
    if[(0=count resultTable)|not 98h=type resultTable; :base];
    okRows:resultTable where (resultTable`status)=`OK;
    stepCount:count resultTable;
    if[0=count okRows; :@[base;`steps;:;stepCount]];
    totalsDict:first 0!select
        totalPnl:sum stepPnl, positionPnlTotal:sum positionPnl, financingTotal:sum financingPnl, txnCostTotal:sum txnCost,
        theoreticalGammaPnlTotal:sum theoreticalGammaPnl, thetaPnlTotal:sum thetaPnl
        from okRows;
    initRow:resultTable 0;
    isCollar:(initRow`mode)=`collar;
    cumPnlSeries:sums okRows`stepPnl;
    maxDrawdownVal:max (maxs cumPnlSeries)-cumPnlSeries;
    spot0:initRow`spot;
    protectionFloor:spot0*1f-stratCfg`putStrikePct;
    capVal:$[isCollar;spot0*1f+stratCfg`callStrikePct;0Nf];
    budgetUsed:$[not isCollar;(initRow`putValue)%spot0;0Nf];
    netPremium:initRow[`putValue]-initRow`callValue;
    totalsDict,`strategyName`mode`steps`netPremium`protectionFloor`cap`budgetUsed`maxDrawdown`status`errorMessage!(
        `collarTailHedge;initRow`mode;stepCount;netPremium;protectionFloor;capVal;budgetUsed;maxDrawdownVal;`OK;"")
 };

.strategy.register[`collarTailHedge;.strategy.collarTailHedge.init;.strategy.collarTailHedge.step;.strategy.collarTailHedge.summary;.strategy.collarTailHedge.defaultConfig];

/ ==================================================================
/ 13. Concrete strategy: putRatioBackspread (1xN)
/ ==================================================================
/ Short 1 near-ATM put + long ratioN OTM puts at a lower strike. Long volatility
/ with downside-skew. Delta-hedge optional via stratCfg.hedgeDelta.

.strategy.putRatioBackspread.__rowEmitCols:`stepIndex`stepDate`spot`shortStrike`longStrike`shortPutPrice`longPutPrice`positionValue`netDelta`hedgePosition`hedgeTrade`txnCost`positionPnl`hedgePnl`financingPnl`thetaPnl`stepPnl`cumulativePnl`theoreticalGammaPnl`status`message;

.strategy.putRatioBackspread.defaultConfig:{[]
    `shortStrikePct`longStrikePct`ratioN`hedgeDelta`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears!(
        0f;0.05;2f;1b;`interval;1;0.05;0f;0f;1f%252f)
 };

.strategy.putRatioBackspread.__validateConfig:{[stratCfg]
    requiredKeys:`shortStrikePct`longStrikePct`ratioN`hedgeDelta`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears;
    missing:requiredKeys where not requiredKeys in key stratCfg;
    if[0<count missing; '"putRatioBackspread config missing keys: ",", " sv string missing];
    if[(stratCfg`longStrikePct)<=stratCfg`shortStrikePct;
        '"putRatioBackspread longStrikePct must exceed shortStrikePct (long strike further OTM)"];
    if[(stratCfg`ratioN)<=1f; '"putRatioBackspread ratioN must be > 1"];
    if[(stratCfg`stepYears)<=0f; '"putRatioBackspread stepYears must be positive"];
 };

.strategy.putRatioBackspread.__priceLeg:{[legTrade;contextDict;model;fdmConfig]
    mktData:.market.createFlatMarketData[contextDict`underlying;contextDict`spot;contextDict`riskFreeRate;contextDict`dividendYield;contextDict`volatility];
    priceRes:.engine.priceOption[legTrade;mktData;model;fdmConfig];
    greeksRes:.greeks.calculateGreeks[legTrade;mktData;model;fdmConfig];
    `unitPrice`unitDelta`unitGamma`unitTheta!(priceRes`unitPrice;first greeksRes`delta;first greeksRes`gamma;first greeksRes`theta)
 };

.strategy.putRatioBackspread.init:{[trade;firstStep;model;fdmConfig;stratCfg]
    .strategy.putRatioBackspread.__validateConfig stratCfg;
    notional:trade`notional;
    spot:firstStep`spot;
    ratioN:stratCfg`ratioN;
    shortStrike:spot*1f-stratCfg`shortStrikePct;
    longStrike:spot*1f-stratCfg`longStrikePct;
    contextDict:`spot`volatility`riskFreeRate`dividendYield`underlying!(
        spot;firstStep`volatility;firstStep`riskFreeRate;firstStep`dividendYield;trade`underlying);
    shortPut:@[trade;(`tradeId;`optionType;`strike);:;(`$(string trade`tradeId),"_SP";`put;shortStrike)];
    longPut:@[trade;(`tradeId;`optionType;`strike);:;(`$(string trade`tradeId),"_LP";`put;longStrike)];
    shortMark:.strategy.putRatioBackspread.__priceLeg[shortPut;contextDict;model;fdmConfig];
    longMark:.strategy.putRatioBackspread.__priceLeg[longPut;contextDict;model;fdmConfig];
    shortUnits:notional;
    longUnits:ratioN*notional;
    shortPrice:shortMark`unitPrice;
    longPrice:longMark`unitPrice;
    positionValue:(longUnits*longPrice)-shortUnits*shortPrice;
    netDelta:(longUnits*longMark`unitDelta)-shortUnits*shortMark`unitDelta;
    netGamma:(longUnits*longMark`unitGamma)-shortUnits*shortMark`unitGamma;
    netTheta:(longUnits*longMark`unitTheta)-shortUnits*shortMark`unitTheta;
    legEntryTxnCost:((longUnits*longPrice)+shortUnits*shortPrice)*stratCfg`txnCostRate;
    hedgeInit:$[stratCfg`hedgeDelta;
        .strategy.__hedgeInit `spot`positionDelta`txnCostRate!(spot;netDelta;stratCfg`txnCostRate);
        `hedgePosition`hedgeTrade`txnCost`cashAdj!(0f;0f;0f;0f)];
    initialTxnCost:legEntryTxnCost+hedgeInit`txnCost;
    initialStepPnl:neg initialTxnCost;
    initialCash:((neg positionValue)-legEntryTxnCost)+hedgeInit`cashAdj;
    rowEmit:.strategy.putRatioBackspread.__rowEmitCols!(
        firstStep`stepIndex;firstStep`stepDate;spot;shortStrike;longStrike;shortPrice;longPrice;
        positionValue;netDelta;hedgeInit`hedgePosition;hedgeInit`hedgeTrade;initialTxnCost;
        0f;0f;0f;0f;initialStepPnl;initialStepPnl;0f;`OK;"");
    `shortPut`longPut`shortUnits`longUnits`shortStrike`longStrike`netPremium`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue`prevPositionGamma`prevPositionTheta`hedgeTrade`txnCost`financingPnl`hedgePnl`positionPnl`thetaPnl`stepPnl`cumulativePnl`rowEmit!(
        shortPut;longPut;shortUnits;longUnits;shortStrike;longStrike;positionValue;initialCash;hedgeInit`hedgePosition;netDelta;1;spot;positionValue;netGamma;netTheta;
        hedgeInit`hedgeTrade;initialTxnCost;0f;0f;0f;0f;initialStepPnl;initialStepPnl;rowEmit)
 };

.strategy.putRatioBackspread.step:{[state;marketStep;trade;model;fdmConfig;stratCfg]
    shortPut:state`shortPut;
    longPut:state`longPut;
    shortUnits:state`shortUnits;
    longUnits:state`longUnits;
    shortStrike:state`shortStrike;
    longStrike:state`longStrike;
    spot:marketStep`spot;
    contextDict:`spot`volatility`riskFreeRate`dividendYield`underlying!(
        spot;marketStep`volatility;marketStep`riskFreeRate;marketStep`dividendYield;trade`underlying);
    shortMark:.strategy.putRatioBackspread.__priceLeg[shortPut;contextDict;model;fdmConfig];
    longMark:.strategy.putRatioBackspread.__priceLeg[longPut;contextDict;model;fdmConfig];
    shortPrice:shortMark`unitPrice;
    longPrice:longMark`unitPrice;
    newPositionValue:(longUnits*longPrice)-shortUnits*shortPrice;
    netDelta:(longUnits*longMark`unitDelta)-shortUnits*shortMark`unitDelta;
    netGamma:(longUnits*longMark`unitGamma)-shortUnits*shortMark`unitGamma;
    netTheta:(longUnits*longMark`unitTheta)-shortUnits*shortMark`unitTheta;
    cashPrev:state`cash;
    financingRate:stratCfg`financingRate;
    hedgeUpdate:$[stratCfg`hedgeDelta;
        .strategy.__hedgeStep[
            `cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue!(
                cashPrev;state`hedgePosition;state`hedgedDelta;state`numRebalances;state`prevSpot;state`prevPositionValue);
            `spot`positionValue`positionDelta`stepIndex`stepYears`txnCostRate`financingRate`rebalanceMode`rebalanceInterval`deltaBand!(
                spot;newPositionValue;netDelta;marketStep`stepIndex;stratCfg`stepYears;stratCfg`txnCostRate;financingRate;stratCfg`rebalanceMode;stratCfg`rebalanceInterval;stratCfg`deltaBand)];
        `cash`hedgePosition`hedgedDelta`numRebalances`hedgeTrade`txnCost`financingPnl`hedgePnl!(
            cashPrev+(financingRate*cashPrev)*stratCfg`stepYears;0f;0f;state`numRebalances;0f;0f;(financingRate*cashPrev)*stratCfg`stepYears;0f)];
    positionPnl:newPositionValue-state`prevPositionValue;
    spotMove:spot-state`prevSpot;
    theoreticalGammaPnl:(0.5*state`prevPositionGamma)*spotMove*spotMove;
    thetaPnl:(state`prevPositionTheta)*stratCfg`stepYears;
    txnCostVal:hedgeUpdate`txnCost;
    stepPnl:(positionPnl+(hedgeUpdate`hedgePnl)+hedgeUpdate`financingPnl)-txnCostVal;
    cumulativePnl:(state`cumulativePnl)+stepPnl;
    rowEmit:.strategy.putRatioBackspread.__rowEmitCols!(
        marketStep`stepIndex;marketStep`stepDate;spot;shortStrike;longStrike;shortPrice;longPrice;
        newPositionValue;netDelta;hedgeUpdate`hedgePosition;hedgeUpdate`hedgeTrade;txnCostVal;
        positionPnl;hedgeUpdate`hedgePnl;hedgeUpdate`financingPnl;thetaPnl;stepPnl;cumulativePnl;theoreticalGammaPnl;`OK;"");
    @[state;`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue`prevPositionGamma`prevPositionTheta`hedgeTrade`txnCost`financingPnl`hedgePnl`positionPnl`thetaPnl`stepPnl`cumulativePnl`rowEmit;:;
        (hedgeUpdate`cash;hedgeUpdate`hedgePosition;hedgeUpdate`hedgedDelta;hedgeUpdate`numRebalances;spot;newPositionValue;netGamma;netTheta;
         hedgeUpdate`hedgeTrade;txnCostVal;hedgeUpdate`financingPnl;hedgeUpdate`hedgePnl;positionPnl;thetaPnl;stepPnl;cumulativePnl;rowEmit)]
 };

.strategy.putRatioBackspread.summary:{[resultTable;stratCfg]
    base:`strategyName`steps`ratioN`shortStrike`longStrike`netPremium`totalPnl`positionPnlTotal`hedgePnlTotal`financingTotal`txnCostTotal`numRebalances`theoreticalGammaPnlTotal`thetaPnlTotal`maxDrawdown`status`errorMessage!(
        `putRatioBackspread;0;stratCfg`ratioN;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0;0Nf;0Nf;0Nf;`ERROR;"empty");
    if[(0=count resultTable)|not 98h=type resultTable; :base];
    okRows:resultTable where (resultTable`status)=`OK;
    stepCount:count resultTable;
    if[0=count okRows; :@[base;`steps;:;stepCount]];
    totalsDict:first 0!select
        totalPnl:sum stepPnl, positionPnlTotal:sum positionPnl, hedgePnlTotal:sum hedgePnl, financingTotal:sum financingPnl, txnCostTotal:sum txnCost,
        theoreticalGammaPnlTotal:sum theoreticalGammaPnl, thetaPnlTotal:sum thetaPnl
        from okRows;
    cumPnl:sums okRows`stepPnl;
    maxDrawdownVal:max (maxs cumPnl)-cumPnl;
    numRebalancesVal:sum 0<>okRows`hedgeTrade;
    initRow:resultTable 0;
    totalsDict,`strategyName`steps`ratioN`shortStrike`longStrike`netPremium`numRebalances`maxDrawdown`status`errorMessage!(
        `putRatioBackspread;stepCount;stratCfg`ratioN;initRow`shortStrike;initRow`longStrike;initRow`positionValue;numRebalancesVal;maxDrawdownVal;`OK;"")
 };

.strategy.register[`putRatioBackspread;.strategy.putRatioBackspread.init;.strategy.putRatioBackspread.step;.strategy.putRatioBackspread.summary;.strategy.putRatioBackspread.defaultConfig];

/ ==================================================================
/ 14. Concrete strategy: ironCondor (defined-risk, 4 legs)
/ ==================================================================
/ Short OTM put (shortPutStrike) + long further-OTM put (longPutStrike) +
/ short OTM call (shortCallStrike) + long further-OTM call (longCallStrike).
/ Held outright (no external delta hedge by default). Net credit collected at
/ entry; analytic max-profit = netCredit; max-loss = widestSpreadWidth - netCredit.
/ Optional IV-rich entry gate: open only when implied vol > entryVolThreshold.

.strategy.ironCondor.__rowEmitCols:`stepIndex`stepDate`spot`shortPutPrice`longPutPrice`shortCallPrice`longCallPrice`positionValue`netDelta`hedgePosition`hedgeTrade`txnCost`positionPnl`hedgePnl`financingPnl`thetaPnl`stepPnl`cumulativePnl`theoreticalGammaPnl`status`message;

.strategy.ironCondor.defaultConfig:{[]
    `shortPutPct`longPutPct`shortCallPct`longCallPct`hedgeDelta`entryVolThreshold`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears!(
        0.05;0.10;0.05;0.10;0b;0f;`interval;1;0.05;0f;0f;1f%252f)
 };

.strategy.ironCondor.__validateConfig:{[stratCfg]
    requiredKeys:`shortPutPct`longPutPct`shortCallPct`longCallPct`hedgeDelta`entryVolThreshold`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears;
    missing:requiredKeys where not requiredKeys in key stratCfg;
    if[0<count missing; '"ironCondor config missing keys: ",", " sv string missing];
    if[(stratCfg`shortPutPct)<=0f; '"ironCondor shortPutPct must be positive"];
    if[(stratCfg`longPutPct)<=stratCfg`shortPutPct;
        '"ironCondor longPutPct must exceed shortPutPct (long put further OTM)"];
    if[(stratCfg`shortCallPct)<=0f; '"ironCondor shortCallPct must be positive"];
    if[(stratCfg`longCallPct)<=stratCfg`shortCallPct;
        '"ironCondor longCallPct must exceed shortCallPct (long call further OTM)"];
    if[(stratCfg`entryVolThreshold)<0f; '"ironCondor entryVolThreshold must be non-negative"];
    if[(stratCfg`stepYears)<=0f; '"ironCondor stepYears must be positive"];
 };

.strategy.ironCondor.__priceLeg:{[legTrade;contextDict;model;fdmConfig]
    mktData:.market.createFlatMarketData[contextDict`underlying;contextDict`spot;contextDict`riskFreeRate;contextDict`dividendYield;contextDict`volatility];
    priceRes:.engine.priceOption[legTrade;mktData;model;fdmConfig];
    greeksRes:.greeks.calculateGreeks[legTrade;mktData;model;fdmConfig];
    `unitPrice`unitDelta`unitGamma`unitTheta!(priceRes`unitPrice;first greeksRes`delta;first greeksRes`gamma;first greeksRes`theta)
 };

.strategy.ironCondor.__flatRow:{[marketStep]
    .strategy.ironCondor.__rowEmitCols!(
        marketStep`stepIndex;marketStep`stepDate;marketStep`spot;
        0Nf;0Nf;0Nf;0Nf;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;
        `flat;"Entry gate closed")
 };

.strategy.ironCondor.init:{[trade;firstStep;model;fdmConfig;stratCfg]
    .strategy.ironCondor.__validateConfig stratCfg;
    notional:trade`notional;
    spot:firstStep`spot;
    vol:firstStep`volatility;
    rfr:firstStep`riskFreeRate;
    divY:firstStep`dividendYield;
    gateOpen:vol>=stratCfg`entryVolThreshold;
    shortPutStrike:spot*1f-stratCfg`shortPutPct;
    longPutStrike:spot*1f-stratCfg`longPutPct;
    shortCallStrike:spot*1f+stratCfg`shortCallPct;
    longCallStrike:spot*1f+stratCfg`longCallPct;
    tidStr:string trade`tradeId;
    shortPutTrade:@[trade;(`tradeId;`optionType;`strike);:;(`$tidStr,"_SP";`put;shortPutStrike)];
    longPutTrade:@[trade;(`tradeId;`optionType;`strike);:;(`$tidStr,"_LP";`put;longPutStrike)];
    shortCallTrade:@[trade;(`tradeId;`optionType;`strike);:;(`$tidStr,"_SC";`call;shortCallStrike)];
    longCallTrade:@[trade;(`tradeId;`optionType;`strike);:;(`$tidStr,"_LC";`call;longCallStrike)];
    putSpreadWidth:shortPutStrike-longPutStrike;
    callSpreadWidth:longCallStrike-shortCallStrike;
    widestSpreadWidth:putSpreadWidth|callSpreadWidth;
    if[not gateOpen;
        flatRow:.strategy.ironCondor.__flatRow firstStep;
        :`gateOpen`notional`shortPutTrade`longPutTrade`shortCallTrade`longCallTrade`shortPutStrike`longPutStrike`shortCallStrike`longCallStrike`putSpreadWidth`callSpreadWidth`widestSpreadWidth`netCreditAtEntry`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue`prevPositionGamma`prevPositionTheta`hedgeTrade`txnCost`financingPnl`hedgePnl`positionPnl`thetaPnl`stepPnl`cumulativePnl`rowEmit!(
            0b;notional;shortPutTrade;longPutTrade;shortCallTrade;longCallTrade;shortPutStrike;longPutStrike;shortCallStrike;longCallStrike;
            putSpreadWidth;callSpreadWidth;widestSpreadWidth;0f;0f;0f;0f;0;spot;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;flatRow)];
    contextDict:`spot`volatility`riskFreeRate`dividendYield`underlying!(
        spot;vol;rfr;divY;trade`underlying);
    spMark:.strategy.ironCondor.__priceLeg[shortPutTrade;contextDict;model;fdmConfig];
    lpMark:.strategy.ironCondor.__priceLeg[longPutTrade;contextDict;model;fdmConfig];
    scMark:.strategy.ironCondor.__priceLeg[shortCallTrade;contextDict;model;fdmConfig];
    lcMark:.strategy.ironCondor.__priceLeg[longCallTrade;contextDict;model;fdmConfig];
    shortPutPrice:spMark`unitPrice;
    longPutPrice:lpMark`unitPrice;
    shortCallPrice:scMark`unitPrice;
    longCallPrice:lcMark`unitPrice;
    netCreditPerUnit:(shortPutPrice+shortCallPrice)-longPutPrice+longCallPrice;
    netCredit:notional*netCreditPerUnit;
    / positionValue: hold +long legs, -short legs (signed mark)
    positionValue:notional*((longPutPrice+longCallPrice)-shortPutPrice+shortCallPrice);
    netDelta:notional*((lpMark`unitDelta)+(lcMark`unitDelta)-(spMark`unitDelta)+scMark`unitDelta);
    netGamma:notional*((lpMark`unitGamma)+(lcMark`unitGamma)-(spMark`unitGamma)+scMark`unitGamma);
    netTheta:notional*((lpMark`unitTheta)+(lcMark`unitTheta)-(spMark`unitTheta)+scMark`unitTheta);
    grossPremium:notional*((shortPutPrice+shortCallPrice)+longPutPrice+longCallPrice);
    legEntryTxnCost:grossPremium*stratCfg`txnCostRate;
    hedgeInit:$[stratCfg`hedgeDelta;
        .strategy.__hedgeInit `spot`positionDelta`txnCostRate!(spot;netDelta;stratCfg`txnCostRate);
        `hedgePosition`hedgeTrade`txnCost`cashAdj!(0f;0f;0f;0f)];
    initialTxnCost:legEntryTxnCost+hedgeInit`txnCost;
    initialStepPnl:neg initialTxnCost;
    initialCash:((neg positionValue)-legEntryTxnCost)+hedgeInit`cashAdj;
    rowEmit:.strategy.ironCondor.__rowEmitCols!(
        firstStep`stepIndex;firstStep`stepDate;spot;shortPutPrice;longPutPrice;shortCallPrice;longCallPrice;
        positionValue;netDelta;hedgeInit`hedgePosition;hedgeInit`hedgeTrade;initialTxnCost;
        0f;0f;0f;0f;initialStepPnl;initialStepPnl;0f;`OK;"");
    `gateOpen`notional`shortPutTrade`longPutTrade`shortCallTrade`longCallTrade`shortPutStrike`longPutStrike`shortCallStrike`longCallStrike`putSpreadWidth`callSpreadWidth`widestSpreadWidth`netCreditAtEntry`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue`prevPositionGamma`prevPositionTheta`hedgeTrade`txnCost`financingPnl`hedgePnl`positionPnl`thetaPnl`stepPnl`cumulativePnl`rowEmit!(
        1b;notional;shortPutTrade;longPutTrade;shortCallTrade;longCallTrade;shortPutStrike;longPutStrike;shortCallStrike;longCallStrike;
        putSpreadWidth;callSpreadWidth;widestSpreadWidth;netCredit;initialCash;hedgeInit`hedgePosition;netDelta;1;spot;positionValue;netGamma;netTheta;
        hedgeInit`hedgeTrade;initialTxnCost;0f;0f;0f;0f;initialStepPnl;initialStepPnl;rowEmit)
 };

.strategy.ironCondor.step:{[state;marketStep;trade;model;fdmConfig;stratCfg]
    if[not state`gateOpen;
        :@[state;`rowEmit;:;.strategy.ironCondor.__flatRow marketStep]];
    notional:state`notional;
    shortPutTrade:state`shortPutTrade;
    longPutTrade:state`longPutTrade;
    shortCallTrade:state`shortCallTrade;
    longCallTrade:state`longCallTrade;
    spot:marketStep`spot;
    contextDict:`spot`volatility`riskFreeRate`dividendYield`underlying!(
        spot;marketStep`volatility;marketStep`riskFreeRate;marketStep`dividendYield;trade`underlying);
    spMark:.strategy.ironCondor.__priceLeg[shortPutTrade;contextDict;model;fdmConfig];
    lpMark:.strategy.ironCondor.__priceLeg[longPutTrade;contextDict;model;fdmConfig];
    scMark:.strategy.ironCondor.__priceLeg[shortCallTrade;contextDict;model;fdmConfig];
    lcMark:.strategy.ironCondor.__priceLeg[longCallTrade;contextDict;model;fdmConfig];
    shortPutPrice:spMark`unitPrice;
    longPutPrice:lpMark`unitPrice;
    shortCallPrice:scMark`unitPrice;
    longCallPrice:lcMark`unitPrice;
    newPositionValue:notional*((longPutPrice+longCallPrice)-shortPutPrice+shortCallPrice);
    netDelta:notional*((lpMark`unitDelta)+(lcMark`unitDelta)-(spMark`unitDelta)+scMark`unitDelta);
    netGamma:notional*((lpMark`unitGamma)+(lcMark`unitGamma)-(spMark`unitGamma)+scMark`unitGamma);
    netTheta:notional*((lpMark`unitTheta)+(lcMark`unitTheta)-(spMark`unitTheta)+scMark`unitTheta);
    cashPrev:state`cash;
    financingRate:stratCfg`financingRate;
    hedgeUpdate:$[stratCfg`hedgeDelta;
        .strategy.__hedgeStep[
            `cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue!(
                cashPrev;state`hedgePosition;state`hedgedDelta;state`numRebalances;state`prevSpot;state`prevPositionValue);
            `spot`positionValue`positionDelta`stepIndex`stepYears`txnCostRate`financingRate`rebalanceMode`rebalanceInterval`deltaBand!(
                spot;newPositionValue;netDelta;marketStep`stepIndex;stratCfg`stepYears;stratCfg`txnCostRate;financingRate;stratCfg`rebalanceMode;stratCfg`rebalanceInterval;stratCfg`deltaBand)];
        `cash`hedgePosition`hedgedDelta`numRebalances`hedgeTrade`txnCost`financingPnl`hedgePnl!(
            cashPrev+(financingRate*cashPrev)*stratCfg`stepYears;0f;0f;state`numRebalances;0f;0f;(financingRate*cashPrev)*stratCfg`stepYears;0f)];
    positionPnl:newPositionValue-state`prevPositionValue;
    spotMove:spot-state`prevSpot;
    theoreticalGammaPnl:(0.5*state`prevPositionGamma)*spotMove*spotMove;
    thetaPnl:(state`prevPositionTheta)*stratCfg`stepYears;
    txnCostVal:hedgeUpdate`txnCost;
    stepPnl:(positionPnl+(hedgeUpdate`hedgePnl)+hedgeUpdate`financingPnl)-txnCostVal;
    cumulativePnl:(state`cumulativePnl)+stepPnl;
    rowEmit:.strategy.ironCondor.__rowEmitCols!(
        marketStep`stepIndex;marketStep`stepDate;spot;shortPutPrice;longPutPrice;shortCallPrice;longCallPrice;
        newPositionValue;netDelta;hedgeUpdate`hedgePosition;hedgeUpdate`hedgeTrade;txnCostVal;
        positionPnl;hedgeUpdate`hedgePnl;hedgeUpdate`financingPnl;thetaPnl;stepPnl;cumulativePnl;theoreticalGammaPnl;`OK;"");
    @[state;`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue`prevPositionGamma`prevPositionTheta`hedgeTrade`txnCost`financingPnl`hedgePnl`positionPnl`thetaPnl`stepPnl`cumulativePnl`rowEmit;:;
        (hedgeUpdate`cash;hedgeUpdate`hedgePosition;hedgeUpdate`hedgedDelta;hedgeUpdate`numRebalances;spot;newPositionValue;netGamma;netTheta;
         hedgeUpdate`hedgeTrade;txnCostVal;hedgeUpdate`financingPnl;hedgeUpdate`hedgePnl;positionPnl;thetaPnl;stepPnl;cumulativePnl;rowEmit)]
 };

.strategy.ironCondor.summary:{[resultTable;stratCfg]
    base:`strategyName`gateOpen`steps`netCredit`maxProfit`maxLoss`lowerBreakeven`upperBreakeven`putSpreadWidth`callSpreadWidth`widestSpreadWidth`totalPnl`positionPnlTotal`hedgePnlTotal`financingTotal`txnCostTotal`theoreticalGammaPnlTotal`thetaPnlTotal`maxDrawdown`status`errorMessage!(
        `ironCondor;0b;0;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;`ERROR;"empty");
    if[(0=count resultTable)|not 98h=type resultTable; :base];
    statusCol:resultTable`status;
    stepCount:count resultTable;
    if[all statusCol=`flat;
        :@[base;(`gateOpen;`steps;`status;`errorMessage;`totalPnl;`positionPnlTotal;`hedgePnlTotal;`financingTotal;`txnCostTotal;`theoreticalGammaPnlTotal;`thetaPnlTotal;`maxDrawdown);:;
            (0b;stepCount;`flat;"Gate closed";0f;0f;0f;0f;0f;0f;0f;0f)]];
    okRows:resultTable where statusCol=`OK;
    if[0=count okRows; :@[base;`steps;:;stepCount]];
    totalsDict:first 0!select
        totalPnl:sum stepPnl, positionPnlTotal:sum positionPnl, hedgePnlTotal:sum hedgePnl, financingTotal:sum financingPnl, txnCostTotal:sum txnCost,
        theoreticalGammaPnlTotal:sum theoreticalGammaPnl, thetaPnlTotal:sum thetaPnl
        from okRows;
    initRow:okRows 0;
    notional:stratCfg`shortPutPct;
    netCreditUnit:((initRow`shortPutPrice)+initRow`shortCallPrice)-(initRow`longPutPrice)+initRow`longCallPrice;
    spot0:initRow`spot;
    shortPutStrike:spot0*1f-stratCfg`shortPutPct;
    longPutStrike:spot0*1f-stratCfg`longPutPct;
    shortCallStrike:spot0*1f+stratCfg`shortCallPct;
    longCallStrike:spot0*1f+stratCfg`longCallPct;
    putSpreadWidth:shortPutStrike-longPutStrike;
    callSpreadWidth:longCallStrike-shortCallStrike;
    widestSpreadWidth:putSpreadWidth|callSpreadWidth;
    netCredit:netCreditUnit;
    maxProfit:netCredit;
    maxLoss:widestSpreadWidth-netCredit;
    lowerBreakeven:shortPutStrike-netCredit;
    upperBreakeven:shortCallStrike+netCredit;
    cumPnl:sums okRows`stepPnl;
    maxDrawdownVal:max (maxs cumPnl)-cumPnl;
    totalsDict,`strategyName`gateOpen`steps`netCredit`maxProfit`maxLoss`lowerBreakeven`upperBreakeven`putSpreadWidth`callSpreadWidth`widestSpreadWidth`maxDrawdown`status`errorMessage!(
        `ironCondor;1b;stepCount;netCredit;maxProfit;maxLoss;lowerBreakeven;upperBreakeven;putSpreadWidth;callSpreadWidth;widestSpreadWidth;maxDrawdownVal;`OK;"")
 };

.strategy.register[`ironCondor;.strategy.ironCondor.init;.strategy.ironCondor.step;.strategy.ironCondor.summary;.strategy.ironCondor.defaultConfig];

/ ==================================================================
/ 15. Concrete strategy: barrierHedge (knock-out + delta-hedge)
/ ==================================================================
/ Hold a long knock-out option (upAndOut call or downAndOut put per trade
/ barrierType) and delta-hedge it on the path. On knock-out event (spot
/ crosses the barrier) the option is settled at the rebate (typically 0),
/ the hedge is closed, and the position is held flat in cash for any
/ remaining steps. Knock detection: upAndOut -> spot >= barrierLevel;
/ downAndOut -> spot <= barrierLevel.

.strategy.barrierHedge.__rowEmitCols:`stepIndex`stepDate`spot`barrierLevel`distanceToBarrier`knockedOut`optionPrice`delta`positionValue`netDelta`hedgePosition`hedgeTrade`txnCost`positionPnl`hedgePnl`financingPnl`thetaPnl`stepPnl`cumulativePnl`theoreticalGammaPnl`status`message;

.strategy.barrierHedge.defaultConfig:{[]
    `rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears!(
        `interval;1;0.05;0f;0f;1f%252f)
 };

.strategy.barrierHedge.__validateConfig:{[stratCfg]
    requiredKeys:`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears;
    missing:requiredKeys where not requiredKeys in key stratCfg;
    if[0<count missing; '"barrierHedge config missing keys: ",", " sv string missing];
    if[not stratCfg[`rebalanceMode] in `interval`band;
        '"barrierHedge rebalanceMode must be interval or band"];
    if[(stratCfg`stepYears)<=0f; '"barrierHedge stepYears must be positive"];
 };

.strategy.barrierHedge.__validateTrade:{[trade]
    if[not `barrierType in key trade; '"barrierHedge requires trade barrierType"];
    if[not (trade`barrierType) in `upAndOut`downAndOut;
        '"barrierHedge requires barrierType upAndOut or downAndOut"];
    if[not `barrierLevel in key trade; '"barrierHedge requires trade barrierLevel"];
    if[(trade`barrierLevel)<=0f; '"barrierHedge barrierLevel must be positive"];
 };

.strategy.barrierHedge.__isKnocked:{[barrierType;spot;barrierLevel]
    $[barrierType=`upAndOut; spot>=barrierLevel; spot<=barrierLevel]
 };

.strategy.barrierHedge.__distanceToBarrier:{[barrierType;spot;barrierLevel]
    $[barrierType=`upAndOut; barrierLevel-spot; spot-barrierLevel]
 };

.strategy.barrierHedge.init:{[trade;firstStep;model;fdmConfig;stratCfg]
    .strategy.barrierHedge.__validateConfig stratCfg;
    .strategy.barrierHedge.__validateTrade trade;
    notional:trade`notional;
    spot:firstStep`spot;
    vol:firstStep`volatility;
    rfr:firstStep`riskFreeRate;
    divY:firstStep`dividendYield;
    barrierType:trade`barrierType;
    barrierLevel:trade`barrierLevel;
    rebate:$[`rebate in key trade; trade`rebate; 0f];
    knockedAtEntry:.strategy.barrierHedge.__isKnocked[barrierType;spot;barrierLevel];
    distanceAtEntry:.strategy.barrierHedge.__distanceToBarrier[barrierType;spot;barrierLevel];
    if[knockedAtEntry;
        flatRow:.strategy.barrierHedge.__rowEmitCols!(
            firstStep`stepIndex;firstStep`stepDate;spot;barrierLevel;distanceAtEntry;1b;
            rebate;0f;notional*rebate;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;
            `flat;"Knocked out at entry");
        :`notional`barrierType`barrierLevel`rebate`knockedOut`knockStep`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue`prevPositionGamma`prevPositionTheta`hedgeTrade`txnCost`financingPnl`hedgePnl`positionPnl`thetaPnl`stepPnl`cumulativePnl`maxAbsHedgeTrade`rowEmit!(
            notional;barrierType;barrierLevel;rebate;1b;firstStep`stepIndex;0f;0f;0f;0;spot;notional*rebate;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;flatRow)];
    mktData:.market.createFlatMarketData[trade`underlying;spot;rfr;divY;vol];
    priceRes:.engine.priceOption[trade;mktData;model;fdmConfig];
    optionPrice:priceRes`unitPrice;
    greeksRes:.greeks.calculateGreeks[trade;mktData;model;fdmConfig];
    deltaVal:first greeksRes`delta;
    gammaVal:first greeksRes`gamma;
    thetaVal:first greeksRes`theta;
    positionValue:notional*optionPrice;
    positionDelta:notional*deltaVal;
    positionGamma:notional*gammaVal;
    positionTheta:notional*thetaVal;
    hedgeInit:.strategy.__hedgeInit `spot`positionDelta`txnCostRate!(spot;positionDelta;stratCfg`txnCostRate);
    initialStepPnl:neg hedgeInit`txnCost;
    initialCash:(neg positionValue)+hedgeInit`cashAdj;
    absHedgeTradeAtEntry:abs hedgeInit`hedgeTrade;
    rowEmit:.strategy.barrierHedge.__rowEmitCols!(
        firstStep`stepIndex;firstStep`stepDate;spot;barrierLevel;distanceAtEntry;0b;
        optionPrice;deltaVal;positionValue;positionDelta;
        hedgeInit`hedgePosition;hedgeInit`hedgeTrade;hedgeInit`txnCost;
        0f;0f;0f;0f;initialStepPnl;initialStepPnl;0f;`OK;"");
    `notional`barrierType`barrierLevel`rebate`knockedOut`knockStep`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue`prevPositionGamma`prevPositionTheta`hedgeTrade`txnCost`financingPnl`hedgePnl`positionPnl`thetaPnl`stepPnl`cumulativePnl`maxAbsHedgeTrade`rowEmit!(
        notional;barrierType;barrierLevel;rebate;0b;0Ni;initialCash;hedgeInit`hedgePosition;positionDelta;1;spot;positionValue;positionGamma;positionTheta;
        hedgeInit`hedgeTrade;hedgeInit`txnCost;0f;0f;0f;0f;initialStepPnl;initialStepPnl;absHedgeTradeAtEntry;rowEmit)
 };

.strategy.barrierHedge.step:{[state;marketStep;trade;model;fdmConfig;stratCfg]
    notional:state`notional;
    barrierType:state`barrierType;
    barrierLevel:state`barrierLevel;
    rebate:state`rebate;
    spot:marketStep`spot;
    distanceVal:.strategy.barrierHedge.__distanceToBarrier[barrierType;spot;barrierLevel];
    cashPrev:state`cash;
    financingRate:stratCfg`financingRate;
    stepYears:stratCfg`stepYears;
    / Already knocked out before this step: only financing flows; everything else flat
    if[state`knockedOut;
        financingPnl:(financingRate*cashPrev)*stepYears;
        newCash:cashPrev+financingPnl;
        stepPnl:financingPnl;
        cumulativePnl:(state`cumulativePnl)+stepPnl;
        rowEmit:.strategy.barrierHedge.__rowEmitCols!(
            marketStep`stepIndex;marketStep`stepDate;spot;barrierLevel;distanceVal;1b;
            rebate;0f;notional*rebate;0f;0f;0f;0f;0f;0f;financingPnl;0f;stepPnl;cumulativePnl;0f;
            `OK;"post knock-out");
        :@[state;`cash`prevSpot`prevPositionValue`hedgePosition`prevPositionGamma`prevPositionTheta`hedgeTrade`txnCost`financingPnl`hedgePnl`positionPnl`thetaPnl`stepPnl`cumulativePnl`rowEmit;:;
            (newCash;spot;notional*rebate;0f;0f;0f;0f;0f;financingPnl;0f;0f;0f;stepPnl;cumulativePnl;rowEmit)]];
    / Detect knock at this step end
    knockNow:.strategy.barrierHedge.__isKnocked[barrierType;spot;barrierLevel];
    if[knockNow;
        / Settle option at rebate; close hedge.
        newPositionValue:notional*rebate;
        prevHedgePos:state`hedgePosition;
        hedgeTrade:neg prevHedgePos;
        txnCostVal:(abs hedgeTrade)*spot*stratCfg`txnCostRate;
        financingPnl:(financingRate*cashPrev)*stepYears;
        positionPnl:newPositionValue-state`prevPositionValue;
        hedgePnl:prevHedgePos*(spot-state`prevSpot);
        stepPnl:(positionPnl+hedgePnl+financingPnl)-txnCostVal;
        newCash:((cashPrev+financingPnl)-hedgeTrade*spot)-txnCostVal;
        cumulativePnl:(state`cumulativePnl)+stepPnl;
        spotMove:spot-state`prevSpot;
        theoreticalGammaPnl:(0.5*state`prevPositionGamma)*spotMove*spotMove;
        thetaPnl:(state`prevPositionTheta)*stepYears;
        maxAbsHedge:(state`maxAbsHedgeTrade)|abs hedgeTrade;
        rowEmit:.strategy.barrierHedge.__rowEmitCols!(
            marketStep`stepIndex;marketStep`stepDate;spot;barrierLevel;distanceVal;1b;
            rebate;0f;newPositionValue;0f;
            0f;hedgeTrade;txnCostVal;
            positionPnl;hedgePnl;financingPnl;thetaPnl;stepPnl;cumulativePnl;theoreticalGammaPnl;
            `OK;"knock-out event");
        :@[state;`knockedOut`knockStep`cash`hedgePosition`hedgedDelta`prevSpot`prevPositionValue`prevPositionGamma`prevPositionTheta`hedgeTrade`txnCost`financingPnl`hedgePnl`positionPnl`thetaPnl`stepPnl`cumulativePnl`maxAbsHedgeTrade`rowEmit;:;
            (1b;marketStep`stepIndex;newCash;0f;0f;spot;newPositionValue;0f;0f;hedgeTrade;txnCostVal;financingPnl;hedgePnl;positionPnl;thetaPnl;stepPnl;cumulativePnl;maxAbsHedge;rowEmit)]];
    / Normal step: still live, re-price barrier option and rebalance hedge
    mktData:.market.createFlatMarketData[trade`underlying;spot;marketStep`riskFreeRate;marketStep`dividendYield;marketStep`volatility];
    priceRes:.engine.priceOption[trade;mktData;model;fdmConfig];
    greeksRes:.greeks.calculateGreeks[trade;mktData;model;fdmConfig];
    optionPrice:priceRes`unitPrice;
    deltaVal:first greeksRes`delta;
    gammaVal:first greeksRes`gamma;
    thetaVal:first greeksRes`theta;
    newPositionValue:notional*optionPrice;
    newPositionDelta:notional*deltaVal;
    newPositionGamma:notional*gammaVal;
    newPositionTheta:notional*thetaVal;
    hedgeState:`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue#state;
    stepInputs:`spot`positionValue`positionDelta`stepIndex`stepYears`txnCostRate`financingRate`rebalanceMode`rebalanceInterval`deltaBand!(
        spot;newPositionValue;newPositionDelta;marketStep`stepIndex;stepYears;stratCfg`txnCostRate;financingRate;stratCfg`rebalanceMode;stratCfg`rebalanceInterval;stratCfg`deltaBand);
    newHedge:.strategy.__hedgeStep[hedgeState;stepInputs];
    spotMove:spot-state`prevSpot;
    theoreticalGammaPnl:(0.5*state`prevPositionGamma)*spotMove*spotMove;
    thetaPnl:(state`prevPositionTheta)*stepYears;
    cumulativePnl:(state`cumulativePnl)+newHedge`stepPnl;
    maxAbsHedge:(state`maxAbsHedgeTrade)|abs newHedge`hedgeTrade;
    rowEmit:.strategy.barrierHedge.__rowEmitCols!(
        marketStep`stepIndex;marketStep`stepDate;spot;barrierLevel;distanceVal;0b;
        optionPrice;deltaVal;newPositionValue;newPositionDelta;
        newHedge`hedgePosition;newHedge`hedgeTrade;newHedge`txnCost;
        newHedge`positionPnl;newHedge`hedgePnl;newHedge`financingPnl;thetaPnl;
        newHedge`stepPnl;cumulativePnl;theoreticalGammaPnl;`OK;"");
    state,newHedge,`prevPositionGamma`prevPositionTheta`cumulativePnl`maxAbsHedgeTrade`rowEmit!(
        newPositionGamma;newPositionTheta;cumulativePnl;maxAbsHedge;rowEmit)
 };

.strategy.barrierHedge.summary:{[resultTable;stratCfg]
    base:`strategyName`steps`knockedOut`knockStep`barrierLevel`maxAbsHedgeTrade`totalPnl`positionPnlTotal`hedgePnlTotal`financingTotal`txnCostTotal`theoreticalGammaPnlTotal`thetaPnlTotal`numRebalances`maxDrawdown`status`errorMessage!(
        `barrierHedge;0;0b;0Ni;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0;0Nf;`ERROR;"empty");
    if[(0=count resultTable)|not 98h=type resultTable; :base];
    okRows:resultTable where (resultTable`status)=`OK;
    stepCount:count resultTable;
    if[0=count okRows; :@[base;`steps;:;stepCount]];
    knockMask:okRows`knockedOut;
    everKnocked:any knockMask;
    knockStepVal:$[everKnocked; first (okRows`stepIndex) where knockMask; 0Ni];
    totalsDict:first 0!select
        totalPnl:sum stepPnl, positionPnlTotal:sum positionPnl, hedgePnlTotal:sum hedgePnl, financingTotal:sum financingPnl, txnCostTotal:sum txnCost,
        theoreticalGammaPnlTotal:sum theoreticalGammaPnl, thetaPnlTotal:sum thetaPnl
        from okRows;
    barrierLevelVal:first okRows`barrierLevel;
    cumPnl:sums okRows`stepPnl;
    maxDrawdownVal:max (maxs cumPnl)-cumPnl;
    numRebalancesVal:sum 0<>okRows`hedgeTrade;
    maxAbsHedgeTradeVal:max abs okRows`hedgeTrade;
    totalsDict,`strategyName`steps`knockedOut`knockStep`barrierLevel`maxAbsHedgeTrade`numRebalances`maxDrawdown`status`errorMessage!(
        `barrierHedge;stepCount;everKnocked;knockStepVal;barrierLevelVal;maxAbsHedgeTradeVal;numRebalancesVal;maxDrawdownVal;`OK;"")
 };

.strategy.register[`barrierHedge;.strategy.barrierHedge.init;.strategy.barrierHedge.step;.strategy.barrierHedge.summary;.strategy.barrierHedge.defaultConfig];

/ ==================================================================
/ 16. Concrete strategy: jumpPremium (real structural model vs BS)
/ ==================================================================
/ At entry: price the option under a jump model (Merton series or Bates MC per
/ stratCfg.jumpModelName) AND under Black-Scholes/FDM. jumpPremium = jumpPx0 - bsPx0.
/ Gate: |jumpPremium| > threshold (absolute or relative to bsPx0 per cfg).
/ Direction: stratCfg.direction = `auto -> sell jump premium (short position) when
/   jumpPx0 > bsPx0, buy (long position) when jumpPx0 < bsPx0;
/ `forceLong / `forceShort override. After entry the position is marked on the path
/ under BS only (for speed); delta-hedged via __hedgeStep. modelDisagreement (vol-bump
/ proxy) is NOT modified; jumpPremium is the structural counterpart.

.strategy.jumpPremium.__rowEmitCols:`stepIndex`stepDate`spot`volatility`optionPrice`delta`positionValue`netDelta`hedgePosition`hedgeTrade`txnCost`positionPnl`hedgePnl`financingPnl`thetaPnl`stepPnl`cumulativePnl`theoreticalGammaPnl`bsPrice0`jumpModelPrice0`jumpPremium`gateOpen`tradeSide`status`message;

.strategy.jumpPremium.defaultConfig:{[]
    `jumpModelName`jumpModelParams`thresholdType`premiumThreshold`direction`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears!(
        `merton;
        `volatility`jumpIntensity`jumpMean`jumpVolatility`riskFreeRate`dividendYield!(0.20;0.5;-0.05;0.20;0.02;0f);
        `absolute;
        0.5;
        `auto;
        `interval;1;0.05;0f;0f;1f%252f)
 };

.strategy.jumpPremium.__validateConfig:{[stratCfg]
    requiredKeys:`jumpModelName`jumpModelParams`thresholdType`premiumThreshold`direction`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears;
    missing:requiredKeys where not requiredKeys in key stratCfg;
    if[0<count missing; '"jumpPremium config missing keys: ",", " sv string missing];
    if[not stratCfg[`jumpModelName] in `merton`bates;
        '"jumpPremium jumpModelName must be merton or bates"];
    if[not stratCfg[`thresholdType] in `absolute`relative;
        '"jumpPremium thresholdType must be absolute or relative"];
    if[not stratCfg[`direction] in `auto`forceLong`forceShort;
        '"jumpPremium direction must be auto, forceLong, or forceShort"];
    if[(stratCfg`premiumThreshold)<0f;
        '"jumpPremium premiumThreshold must be non-negative"];
    if[(stratCfg`stepYears)<=0f; '"jumpPremium stepYears must be positive"];
 };

.strategy.jumpPremium.__priceJumpModel:{[trade;mktData;stratCfg]
    jumpName:stratCfg`jumpModelName;
    jumpParams:stratCfg`jumpModelParams;
    if[jumpName=`merton;
        configDict:`mertonParams`pricingMethod`termCount!(jumpParams;`series;30);
        priceRes:.merton.priceEuropean[trade;mktData;configDict];
        :priceRes`unitPrice];
    if[jumpName=`bates;
        mcConfig:$[`mcConfig in key stratCfg; stratCfg`mcConfig; .montecarlo.defaultMcConfig[]];
        configDict:`batesParams`mcConfig!(jumpParams;mcConfig);
        priceRes:.bates.priceEuropean[trade;mktData;configDict];
        :priceRes`unitPrice];
    '"jumpPremium unsupported jumpModelName: ",string jumpName
 };

.strategy.jumpPremium.__resolveSide:{[direction;jumpPx;bsPx]
    if[direction=`forceLong; :`long];
    if[direction=`forceShort; :`short];
    $[jumpPx>bsPx; `short; `long]
 };

.strategy.jumpPremium.__flatRow:{[marketStep;tradeSideVal;bsPx0;jumpPx0;premium]
    .strategy.jumpPremium.__rowEmitCols!(
        marketStep`stepIndex;marketStep`stepDate;marketStep`spot;marketStep`volatility;
        0Nf;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;bsPx0;jumpPx0;premium;0b;tradeSideVal;
        `flat;"Premium below threshold")
 };

.strategy.jumpPremium.init:{[trade;firstStep;model;fdmConfig;stratCfg]
    .strategy.jumpPremium.__validateConfig stratCfg;
    notional:trade`notional;
    spot:firstStep`spot;
    vol:firstStep`volatility;
    rfr:firstStep`riskFreeRate;
    divY:firstStep`dividendYield;
    mktData:.market.createFlatMarketData[trade`underlying;spot;rfr;divY;vol];
    bsPx:(.engine.priceOption[trade;mktData;model;fdmConfig])`unitPrice;
    jumpPx:.strategy.jumpPremium.__priceJumpModel[trade;mktData;stratCfg];
    jumpPremium:jumpPx-bsPx;
    threshold:$[stratCfg[`thresholdType]=`absolute;
        stratCfg`premiumThreshold;
        (stratCfg`premiumThreshold)*abs bsPx];
    gateOpen:threshold<abs jumpPremium;
    tradeSideVal:.strategy.jumpPremium.__resolveSide[stratCfg`direction;jumpPx;bsPx];
    if[not gateOpen;
        flatRow:.strategy.jumpPremium.__flatRow[firstStep;tradeSideVal;bsPx;jumpPx;jumpPremium];
        :`gateOpen`notional`tradeSide`bsPrice0`jumpModelPrice0`jumpPremium`jumpModelName`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue`prevPositionGamma`prevPositionTheta`hedgeTrade`txnCost`financingPnl`hedgePnl`positionPnl`thetaPnl`stepPnl`cumulativePnl`rowEmit!(
            0b;notional;tradeSideVal;bsPx;jumpPx;jumpPremium;stratCfg`jumpModelName;0f;0f;0f;0;spot;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;0f;flatRow)];
    sideSign:$[tradeSideVal=`long; 1f; -1f];
    optionUnits:notional*sideSign;
    greeksRes:.greeks.calculateGreeks[trade;mktData;model;fdmConfig];
    deltaVal:first greeksRes`delta;
    gammaVal:first greeksRes`gamma;
    thetaVal:first greeksRes`theta;
    positionValue:optionUnits*bsPx;
    positionDelta:optionUnits*deltaVal;
    positionGamma:optionUnits*gammaVal;
    positionTheta:optionUnits*thetaVal;
    hedgeInit:.strategy.__hedgeInit `spot`positionDelta`txnCostRate!(spot;positionDelta;stratCfg`txnCostRate);
    initialCash:(neg positionValue)+hedgeInit`cashAdj;
    initialStepPnl:neg hedgeInit`txnCost;
    rowEmit:.strategy.jumpPremium.__rowEmitCols!(
        firstStep`stepIndex;firstStep`stepDate;spot;vol;bsPx;deltaVal;positionValue;positionDelta;
        hedgeInit`hedgePosition;hedgeInit`hedgeTrade;hedgeInit`txnCost;
        0f;0f;0f;0f;initialStepPnl;initialStepPnl;0f;bsPx;jumpPx;jumpPremium;1b;tradeSideVal;`OK;"");
    `gateOpen`notional`optionUnits`tradeSide`bsPrice0`jumpModelPrice0`jumpPremium`jumpModelName`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue`prevPositionGamma`prevPositionTheta`hedgeTrade`txnCost`financingPnl`hedgePnl`positionPnl`thetaPnl`stepPnl`cumulativePnl`rowEmit!(
        1b;notional;optionUnits;tradeSideVal;bsPx;jumpPx;jumpPremium;stratCfg`jumpModelName;initialCash;hedgeInit`hedgePosition;positionDelta;1;spot;positionValue;positionGamma;positionTheta;
        hedgeInit`hedgeTrade;hedgeInit`txnCost;0f;0f;0f;0f;initialStepPnl;initialStepPnl;rowEmit)
 };

.strategy.jumpPremium.step:{[state;marketStep;trade;model;fdmConfig;stratCfg]
    if[not state`gateOpen;
        :@[state;`rowEmit;:;.strategy.jumpPremium.__flatRow[marketStep;state`tradeSide;state`bsPrice0;state`jumpModelPrice0;state`jumpPremium]]];
    optionUnits:state`optionUnits;
    spot:marketStep`spot;
    vol:marketStep`volatility;
    mktData:.market.createFlatMarketData[trade`underlying;spot;marketStep`riskFreeRate;marketStep`dividendYield;vol];
    bsPx:(.engine.priceOption[trade;mktData;model;fdmConfig])`unitPrice;
    greeksRes:.greeks.calculateGreeks[trade;mktData;model;fdmConfig];
    deltaVal:first greeksRes`delta;
    gammaVal:first greeksRes`gamma;
    thetaVal:first greeksRes`theta;
    newPositionValue:optionUnits*bsPx;
    newPositionDelta:optionUnits*deltaVal;
    newPositionGamma:optionUnits*gammaVal;
    newPositionTheta:optionUnits*thetaVal;
    hedgeState:`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue#state;
    stepInputs:`spot`positionValue`positionDelta`stepIndex`stepYears`txnCostRate`financingRate`rebalanceMode`rebalanceInterval`deltaBand!(
        spot;newPositionValue;newPositionDelta;marketStep`stepIndex;stratCfg`stepYears;stratCfg`txnCostRate;stratCfg`financingRate;stratCfg`rebalanceMode;stratCfg`rebalanceInterval;stratCfg`deltaBand);
    newHedge:.strategy.__hedgeStep[hedgeState;stepInputs];
    spotMove:spot-state`prevSpot;
    theoreticalGammaPnl:(0.5*state`prevPositionGamma)*spotMove*spotMove;
    thetaPnl:(state`prevPositionTheta)*stratCfg`stepYears;
    cumulativePnl:(state`cumulativePnl)+newHedge`stepPnl;
    rowEmit:.strategy.jumpPremium.__rowEmitCols!(
        marketStep`stepIndex;marketStep`stepDate;spot;vol;bsPx;deltaVal;newPositionValue;newPositionDelta;
        newHedge`hedgePosition;newHedge`hedgeTrade;newHedge`txnCost;
        newHedge`positionPnl;newHedge`hedgePnl;newHedge`financingPnl;thetaPnl;
        newHedge`stepPnl;cumulativePnl;theoreticalGammaPnl;state`bsPrice0;state`jumpModelPrice0;state`jumpPremium;1b;state`tradeSide;`OK;"");
    state,newHedge,`prevPositionGamma`prevPositionTheta`cumulativePnl`rowEmit!(
        newPositionGamma;newPositionTheta;cumulativePnl;rowEmit)
 };

.strategy.jumpPremium.summary:{[resultTable;stratCfg]
    base:`strategyName`gateOpen`tradeSide`jumpModelName`bsPrice0`jumpModelPrice0`jumpPremium`steps`totalPnl`positionPnlTotal`hedgePnlTotal`financingTotal`txnCostTotal`theoreticalGammaPnlTotal`thetaPnlTotal`numRebalances`maxDrawdown`status`errorMessage!(
        `jumpPremium;0b;`auto;stratCfg`jumpModelName;0Nf;0Nf;0Nf;0;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0;0Nf;`ERROR;"empty");
    if[(0=count resultTable)|not 98h=type resultTable; :base];
    statusCol:resultTable`status;
    stepCount:count resultTable;
    initRow:resultTable 0;
    if[all statusCol=`flat;
        :@[base;(`gateOpen;`tradeSide;`bsPrice0;`jumpModelPrice0;`jumpPremium;`steps;`status;`errorMessage;`totalPnl;`positionPnlTotal;`hedgePnlTotal;`financingTotal;`txnCostTotal;`theoreticalGammaPnlTotal;`thetaPnlTotal;`maxDrawdown);:;
            (0b;initRow`tradeSide;initRow`bsPrice0;initRow`jumpModelPrice0;initRow`jumpPremium;stepCount;`flat;"Gate closed";0f;0f;0f;0f;0f;0f;0f;0f)]];
    okRows:resultTable where statusCol=`OK;
    if[0=count okRows; :@[base;`steps;:;stepCount]];
    totalsDict:first 0!select
        totalPnl:sum stepPnl, positionPnlTotal:sum positionPnl, hedgePnlTotal:sum hedgePnl, financingTotal:sum financingPnl, txnCostTotal:sum txnCost,
        theoreticalGammaPnlTotal:sum theoreticalGammaPnl, thetaPnlTotal:sum thetaPnl
        from okRows;
    cumPnl:sums okRows`stepPnl;
    maxDrawdownVal:max (maxs cumPnl)-cumPnl;
    numRebalancesVal:sum 0<>okRows`hedgeTrade;
    initRowOk:okRows 0;
    totalsDict,`strategyName`gateOpen`tradeSide`jumpModelName`bsPrice0`jumpModelPrice0`jumpPremium`steps`numRebalances`maxDrawdown`status`errorMessage!(
        `jumpPremium;1b;initRowOk`tradeSide;stratCfg`jumpModelName;initRowOk`bsPrice0;initRowOk`jumpModelPrice0;initRowOk`jumpPremium;stepCount;numRebalancesVal;maxDrawdownVal;`OK;"")
 };

.strategy.register[`jumpPremium;.strategy.jumpPremium.init;.strategy.jumpPremium.step;.strategy.jumpPremium.summary;.strategy.jumpPremium.defaultConfig];

/ ==================================================================
/ 17. Path adapter: fromCorrelated (multi-asset, correlated GBM)
/ ==================================================================
/ Deterministic Cholesky-correlated GBM paths across n names. Returns a
/ bundle with per-name standard path tables, the derived index path
/ (weighted basket level), weights, vols, and the input correlation matrix.

.strategy.path.fromCorrelated:{[pathCfg]
    required:`names`weights`spot0`vols`correlationMatrix`steps`stepYears`riskFreeRate`dividendYields`seed;
    missing:required where not required in key pathCfg;
    if[0<count missing; '"strategy.path.fromCorrelated missing keys: ",", " sv string missing];
    namesList:pathCfg`names;
    weights:`float$pathCfg`weights;
    spot0Vec:`float$pathCfg`spot0;
    volsVec:`float$pathCfg`vols;
    corrMat:pathCfg`correlationMatrix;
    stepsTotal:pathCfg`steps;
    dt:`float$pathCfg`stepYears;
    rfr:`float$pathCfg`riskFreeRate;
    divYields:`float$pathCfg`dividendYields;
    seed:pathCfg`seed;
    n:count namesList;
    if[stepsTotal<=1; '"strategy.path.fromCorrelated steps must be > 1"];
    if[dt<=0f; '"strategy.path.fromCorrelated stepYears must be positive"];
    if[n<2; '"strategy.path.fromCorrelated need at least 2 names"];
    if[n<>count weights; '"strategy.path.fromCorrelated weights length mismatch"];
    if[n<>count spot0Vec; '"strategy.path.fromCorrelated spot0 length mismatch"];
    if[n<>count volsVec; '"strategy.path.fromCorrelated vols length mismatch"];
    if[n<>count divYields; '"strategy.path.fromCorrelated dividendYields length mismatch"];
    if[n<>count corrMat; '"strategy.path.fromCorrelated correlationMatrix not n x n"];
    if[not all n=count each corrMat; '"strategy.path.fromCorrelated correlationMatrix not square"];
    if[not .correlation.isSymmetric[corrMat;1e-10]; '"strategy.path.fromCorrelated correlationMatrix not symmetric"];
    if[1e-8<abs (sum weights)-1f; '"strategy.path.fromCorrelated weights must sum to 1"];
    if[any volsVec<0f; '"strategy.path.fromCorrelated vols must be non-negative"];
    if[any spot0Vec<=0f; '"strategy.path.fromCorrelated spot0 must be positive"];
    if[any divYields<0f; '"strategy.path.fromCorrelated dividendYields must be non-negative"];
    cholL:.correlation.__cholesky corrMat;
    nIncr:stepsTotal-1;
    flatNormals:.montecarlo.__generateNormals[n*nIncr;seed];
    rawMat:nIncr cut flatNormals;
    correlatedFull:cholL mmu rawMat;
    sqrtDt:sqrt dt;
    halfVarVec:0.5*volsVec*volsVec;
    fullDrift:((rfr-divYields)-halfVarVec)*dt;
    startDate:$[`startDate in key pathCfg; pathCfg`startDate; 2024.01.01];
    stepDates:startDate+til stepsTotal;
    nameIdx:0;
    perNameTables:();
    while[nameIdx<n;
        zN:correlatedFull[nameIdx];
        incr:fullDrift[nameIdx]+volsVec[nameIdx]*sqrtDt*zN;
        logSpotInc:(log spot0Vec nameIdx)+sums incr;
        spots:spot0Vec[nameIdx],exp logSpotInc;
        tbl:flip .strategy.path.__schemaCols!(
            til stepsTotal;
            stepDates;
            spots;
            stepsTotal#volsVec nameIdx;
            stepsTotal#rfr;
            stepsTotal#divYields nameIdx;
            stepsTotal#0Nf;
            stepsTotal#`OK);
        perNameTables,:enlist tbl;
        nameIdx+:1];
    pathTables:namesList!perNameTables;
    spotMatrix:perNameTables[;`spot];
    indexSpots:sum weights*spotMatrix;
    avgDivY:weights mmu divYields;
    indexPath:flip .strategy.path.__schemaCols!(
        til stepsTotal;
        stepDates;
        indexSpots;
        stepsTotal#0Nf;
        stepsTotal#rfr;
        stepsTotal#avgDivY;
        stepsTotal#0Nf;
        stepsTotal#`OK);
    `names`weights`pathTables`indexPath`correlationMatrix`vols!(
        namesList;weights;pathTables;indexPath;corrMat;volsVec)
 };

/ ==================================================================
/ 18. Concrete strategy: dispersion (index vs constituents)
/ ==================================================================
/ Short index ATM straddle + long weighted constituent ATM straddles (or
/ reverse). Holds the index straddle UNHEDGED (correlation exposure). Each
/ constituent straddle is delta-hedged in its own underlying via __hedgeStep.
/ Per-underlying hedge book is a TABLE keyed by name. The pathTable passed to
/ .strategy.run IS the index path; per-name spots/vols come from
/ stratCfg.pathBundle.

.strategy.dispersion.__rowEmitCols:`stepIndex`stepDate`spot`indexStraddleValue`constituentStraddleSum`positionValue`netDelta`hedgeNotional`txnCost`positionPnl`hedgePnl`financingPnl`stepPnl`cumulativePnl`status`message;

.strategy.dispersion.defaultConfig:{[]
    `pathBundle`indexVol`direction`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears!(
        ()!();
        0.20;
        `shortIndex;
        `interval;1;0.05;0f;0f;1f%252f)
 };

.strategy.dispersion.__validateConfig:{[stratCfg]
    requiredKeys:`pathBundle`indexVol`direction`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears;
    missing:requiredKeys where not requiredKeys in key stratCfg;
    if[0<count missing; '"dispersion config missing keys: ",", " sv string missing];
    if[not stratCfg[`direction] in `shortIndex`longIndex;
        '"dispersion direction must be shortIndex or longIndex"];
    if[(stratCfg`indexVol)<=0f; '"dispersion indexVol must be positive"];
    if[(stratCfg`stepYears)<=0f; '"dispersion stepYears must be positive"];
    bundle:stratCfg`pathBundle;
    bundleKeys:`names`weights`pathTables`indexPath`correlationMatrix`vols;
    bundleMissing:bundleKeys where not bundleKeys in key bundle;
    if[0<count bundleMissing; '"dispersion pathBundle missing keys: ",", " sv string bundleMissing];
 };

.strategy.dispersion.__priceStraddle:{[straddleCtx;model;fdmConfig]
    spot:straddleCtx`spot;
    strike:straddleCtx`strike;
    timeToExpiry:straddleCtx`timeToExpiry;
    rfr:straddleCtx`riskFreeRate;
    divY:straddleCtx`dividendYield;
    vol:straddleCtx`volatility;
    underlyingSym:straddleCtx`underlying;
    callTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
        `$"_C_",string underlyingSym;underlyingSym;`equityOption;`european;`call;strike;timeToExpiry;1f);
    putTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
        `$"_P_",string underlyingSym;underlyingSym;`equityOption;`european;`put;strike;timeToExpiry;1f);
    mktData:.market.createFlatMarketData[underlyingSym;spot;rfr;divY;vol];
    callPx:(.engine.priceOption[callTrade;mktData;model;fdmConfig])`unitPrice;
    putPx:(.engine.priceOption[putTrade;mktData;model;fdmConfig])`unitPrice;
    callGreeks:.greeks.calculateGreeks[callTrade;mktData;model;fdmConfig];
    putGreeks:.greeks.calculateGreeks[putTrade;mktData;model;fdmConfig];
    straddleValue:callPx+putPx;
    straddleDelta:(first callGreeks`delta)+first putGreeks`delta;
    `value`delta!(straddleValue;straddleDelta)
 };

.strategy.dispersion.init:{[trade;firstStep;model;fdmConfig;stratCfg]
    .strategy.dispersion.__validateConfig stratCfg;
    bundle:stratCfg`pathBundle;
    namesList:bundle`names;
    weights:bundle`weights;
    nameVols:bundle`vols;
    pathTables:bundle`pathTables;
    n:count namesList;
    notional:trade`notional;
    direction:stratCfg`direction;
    indexSign:$[direction=`shortIndex; -1f; 1f];
    constituentSign:neg indexSign;
    indexVol:stratCfg`indexVol;
    indexSpot:firstStep`spot;
    rfr:firstStep`riskFreeRate;
    indexDivY:firstStep`dividendYield;
    expiryYears:trade`expiry;
    idxCtxInit:`spot`strike`timeToExpiry`riskFreeRate`dividendYield`volatility`underlying!(
        indexSpot;indexSpot;expiryYears;rfr;indexDivY;indexVol;trade`underlying);
    indexStraddleMark:.strategy.dispersion.__priceStraddle[idxCtxInit;model;fdmConfig];
    perNameMarks:();
    perNameSpots0:();
    perNameDivYields:();
    nameIdx:0;
    while[nameIdx<n;
        nm:namesList nameIdx;
        nmTbl:pathTables nm;
        nmRow0:nmTbl 0;
        nmSpot0:nmRow0`spot;
        nmDivY:nmRow0`dividendYield;
        nmVol:nameVols nameIdx;
        nmCtx:`spot`strike`timeToExpiry`riskFreeRate`dividendYield`volatility`underlying!(
            nmSpot0;nmSpot0;expiryYears;rfr;nmDivY;nmVol;nm);
        mk:.strategy.dispersion.__priceStraddle[nmCtx;model;fdmConfig];
        perNameMarks,:enlist mk;
        perNameSpots0,:nmSpot0;
        perNameDivYields,:nmDivY;
        nameIdx+:1];
    constituentValues:perNameMarks[;`value];
    constituentDeltas:perNameMarks[;`delta];
    constituentStraddleSum:notional*sum weights*constituentValues;
    indexStraddleValue:notional*indexStraddleMark`value;
    positionValue:(indexSign*indexStraddleValue)+constituentSign*constituentStraddleSum;
    netConstituentDeltas:notional*constituentSign*weights*constituentDeltas;
    indexHedgeNetDeltaContribution:notional*indexSign*indexStraddleMark`delta;
    netDelta:indexHedgeNetDeltaContribution+sum netConstituentDeltas;
    / Per-name hedge book (table). Hedge each constituent leg's delta.
    hedgePositions:neg netConstituentDeltas;
    hedgeTradesInit:hedgePositions;
    txnCostInit:sum (abs hedgeTradesInit)*perNameSpots0*stratCfg`txnCostRate;
    legEntryCash:neg positionValue;
    hedgeEntryCash:neg (sum hedgePositions*perNameSpots0)+txnCostInit;
    initialCash:legEntryCash+hedgeEntryCash;
    initialStepPnl:neg txnCostInit;
    hedgeBook:flip `name`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue!(
        namesList;
        n#0f;
        hedgePositions;
        constituentSign*weights*notional*constituentDeltas;
        n#1;
        perNameSpots0;
        notional*constituentSign*weights*constituentValues);
    rowEmit:.strategy.dispersion.__rowEmitCols!(
        firstStep`stepIndex;firstStep`stepDate;indexSpot;
        notional*indexStraddleMark`value;constituentStraddleSum;positionValue;netDelta;
        sum abs hedgePositions*perNameSpots0;txnCostInit;
        0f;0f;0f;initialStepPnl;initialStepPnl;`OK;"");
    `notional`direction`indexSign`constituentSign`indexVol`expiryYears`pathBundle`hedgeBook`cash`prevSpot`prevPositionValue`indexStraddleMarkPrev`constituentStraddleSumPrev`prevSpots`rowEmit!(
        notional;direction;indexSign;constituentSign;indexVol;expiryYears;bundle;hedgeBook;initialCash;indexSpot;positionValue;indexStraddleMark`value;constituentStraddleSum;perNameSpots0;rowEmit)
 };

.strategy.dispersion.step:{[state;marketStep;trade;model;fdmConfig;stratCfg]
    bundle:state`pathBundle;
    namesList:bundle`names;
    weights:bundle`weights;
    nameVols:bundle`vols;
    pathTables:bundle`pathTables;
    n:count namesList;
    notional:state`notional;
    indexSign:state`indexSign;
    constituentSign:state`constituentSign;
    indexVol:state`indexVol;
    expiryYears:state`expiryYears;
    indexSpot:marketStep`spot;
    rfr:marketStep`riskFreeRate;
    indexDivY:marketStep`dividendYield;
    stepIndex:marketStep`stepIndex;
    / Strike was set at entry to entry-time index spot.
    indexStrike0:first (bundle`indexPath)[;`spot];
    idxCtx:`spot`strike`timeToExpiry`riskFreeRate`dividendYield`volatility`underlying!(
        indexSpot;indexStrike0;expiryYears;rfr;indexDivY;indexVol;trade`underlying);
    indexStraddleMark:.strategy.dispersion.__priceStraddle[idxCtx;model;fdmConfig];
    perNameMarks:();
    perNameSpots:();
    nameIdx:0;
    while[nameIdx<n;
        nm:namesList nameIdx;
        nmTbl:pathTables nm;
        nmRow:nmTbl stepIndex;
        nmSpot:nmRow`spot;
        nmDivY:nmRow`dividendYield;
        nmVol:nameVols nameIdx;
        nmStrike0:first (nmTbl)[;`spot];
        nmCtx:`spot`strike`timeToExpiry`riskFreeRate`dividendYield`volatility`underlying!(
            nmSpot;nmStrike0;expiryYears;rfr;nmDivY;nmVol;nm);
        mk:.strategy.dispersion.__priceStraddle[nmCtx;model;fdmConfig];
        perNameMarks,:enlist mk;
        perNameSpots,:nmSpot;
        nameIdx+:1];
    constituentValues:perNameMarks[;`value];
    constituentDeltas:perNameMarks[;`delta];
    constituentStraddleSum:notional*sum weights*constituentValues;
    indexStraddleValue:notional*indexStraddleMark`value;
    newPositionValue:(indexSign*indexStraddleValue)+constituentSign*constituentStraddleSum;
    netConstituentDeltas:notional*constituentSign*weights*constituentDeltas;
    indexHedgeNetDeltaContribution:notional*indexSign*indexStraddleMark`delta;
    netDelta:indexHedgeNetDeltaContribution+sum netConstituentDeltas;
    hedgeBookPrev:state`hedgeBook;
    cashPrev:state`cash;
    financingPnl:(stratCfg`financingRate)*cashPrev*stratCfg`stepYears;
    / For each name, compute hedge step (treating netConstituentDelta as positionDelta)
    txnCostRate:stratCfg`txnCostRate;
    rebalanceMode:stratCfg`rebalanceMode;
    rebalanceInterval:stratCfg`rebalanceInterval;
    deltaBand:stratCfg`deltaBand;
    nameIdx:0;
    hedgeRowsNew:();
    aggPositionPnl:0f;
    aggHedgePnl:0f;
    aggTxnCost:0f;
    aggHedgeNotional:0f;
    while[nameIdx<n;
        nm:namesList nameIdx;
        prevRow:first 0!hedgeBookPrev where (hedgeBookPrev`name)=nm;
        legPositionValue:notional*constituentSign*weights[nameIdx]*constituentValues nameIdx;
        legPositionDelta:netConstituentDeltas nameIdx;
        nmSpot:perNameSpots nameIdx;
        prevNmSpot:state[`prevSpots] nameIdx;
        positionPnlLeg:legPositionValue-prevRow`prevPositionValue;
        hedgePnlLeg:(prevRow`hedgePosition)*nmSpot-prevNmSpot;
        intervalTrigger:0=stepIndex mod rebalanceInterval;
        bandTrigger:(abs legPositionDelta-prevRow`hedgedDelta)>deltaBand;
        shouldRebalance:$[rebalanceMode=`interval; intervalTrigger; bandTrigger];
        newHedgePos:$[shouldRebalance; neg legPositionDelta; prevRow`hedgePosition];
        hedgeTrade:newHedgePos-prevRow`hedgePosition;
        txnCostLeg:(abs hedgeTrade)*nmSpot*txnCostRate;
        newHedgedDelta:$[shouldRebalance; legPositionDelta; prevRow`hedgedDelta];
        rebalanceIncrement:$[shouldRebalance&0<>hedgeTrade; 1; 0];
        numReb:(prevRow`numRebalances)+rebalanceIncrement;
        / per-leg cash flow: hedge trade cost + txn cost
        legHedgeCash:neg ((hedgeTrade*nmSpot)+txnCostLeg);
        newLegCash:(prevRow`cash)+legHedgeCash;
        hedgeRowsNew,:enlist `name`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue!(
            nm;newLegCash;newHedgePos;newHedgedDelta;numReb;nmSpot;legPositionValue);
        aggPositionPnl+:positionPnlLeg;
        aggHedgePnl+:hedgePnlLeg;
        aggTxnCost+:txnCostLeg;
        aggHedgeNotional+:abs newHedgePos*nmSpot;
        nameIdx+:1];
    newHedgeBook:hedgeRowsNew;
    / Index straddle position P&L
    indexPositionPnlLeg:(indexSign*indexStraddleValue)-(state`indexStraddleMarkPrev)*indexSign*notional;
    aggPositionPnl+:indexPositionPnlLeg;
    / Total cash update: financing on cashPrev plus all leg hedge cash flows (carried in hedgeBook now)
    legCashChange:sum (newHedgeBook`cash)-hedgeBookPrev`cash;
    newCash:(cashPrev+financingPnl)+legCashChange;
    stepPnl:(aggPositionPnl+aggHedgePnl+financingPnl)-aggTxnCost;
    cumulativePnl:(state[`rowEmit]`cumulativePnl)+stepPnl;
    rowEmit:.strategy.dispersion.__rowEmitCols!(
        stepIndex;marketStep`stepDate;indexSpot;
        indexStraddleValue;constituentStraddleSum;newPositionValue;netDelta;
        aggHedgeNotional;aggTxnCost;
        aggPositionPnl;aggHedgePnl;financingPnl;stepPnl;cumulativePnl;`OK;"");
    @[state;`hedgeBook`cash`prevSpot`prevPositionValue`indexStraddleMarkPrev`constituentStraddleSumPrev`prevSpots`rowEmit;:;
        (newHedgeBook;newCash;indexSpot;newPositionValue;indexStraddleMark`value;constituentStraddleSum;perNameSpots;rowEmit)]
 };

.strategy.dispersion.summary:{[resultTable;stratCfg]
    base:`strategyName`steps`direction`numConstituents`indexVol`impliedCorrAtEntry`realizedCorr`correlationPremium`totalPnl`positionPnlTotal`hedgePnlTotal`financingTotal`txnCostTotal`maxDrawdown`status`errorMessage!(
        `dispersion;0;stratCfg`direction;0;stratCfg`indexVol;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;`ERROR;"empty");
    if[(0=count resultTable)|not 98h=type resultTable; :base];
    okRows:resultTable where (resultTable`status)=`OK;
    stepCount:count resultTable;
    if[0=count okRows; :@[base;`steps;:;stepCount]];
    totalsDict:first 0!select
        totalPnl:sum stepPnl, positionPnlTotal:sum positionPnl, hedgePnlTotal:sum hedgePnl, financingTotal:sum financingPnl, txnCostTotal:sum txnCost
        from okRows;
    bundle:stratCfg`pathBundle;
    namesList:bundle`names;
    weights:bundle`weights;
    nameVols:bundle`vols;
    corrMat:bundle`correlationMatrix;
    pathTables:bundle`pathTables;
    n:count namesList;
    indexVolUsed:stratCfg`indexVol;
    sumWiSqSigSq:sum (weights*weights)*nameVols*nameVols;
    crossTerm:0f;
    iIdx:0;
    while[iIdx<n;
        jIdx:iIdx+1;
        while[jIdx<n;
            crossTerm+:weights[iIdx]*weights[jIdx]*nameVols[iIdx]*nameVols[jIdx];
            jIdx+:1];
        iIdx+:1];
    impliedAvgCorr:$[crossTerm>0f;
        ((indexVolUsed*indexVolUsed)-sumWiSqSigSq)%(2f*crossTerm);
        0Nf];
    / Realized: pairwise log-return correlations averaged over off-diagonal pairs
    pathSpotsList:();
    pathReturnsList:();
    nameIdxSp:0;
    while[nameIdxSp<n;
        sp:(pathTables namesList[nameIdxSp])`spot;
        logS:log sp;
        rets:1_logS-prev logS;
        pathSpotsList,:enlist sp;
        pathReturnsList,:enlist rets;
        nameIdxSp+:1];
    pairwiseCorrs:();
    iIdx:0;
    while[iIdx<n;
        jIdx:iIdx+1;
        while[jIdx<n;
            xV:pathReturnsList iIdx;
            yV:pathReturnsList jIdx;
            pairwiseCorrs,:enlist xV cor yV;
            jIdx+:1];
        iIdx+:1];
    realizedAvgCorr:$[0<count pairwiseCorrs; avg pairwiseCorrs; 0Nf];
    correlationPremium:impliedAvgCorr-realizedAvgCorr;
    cumPnl:sums okRows`stepPnl;
    maxDrawdownVal:max (maxs cumPnl)-cumPnl;
    totalsDict,`strategyName`steps`direction`numConstituents`indexVol`impliedCorrAtEntry`realizedCorr`correlationPremium`maxDrawdown`status`errorMessage!(
        `dispersion;stepCount;stratCfg`direction;n;indexVolUsed;impliedAvgCorr;realizedAvgCorr;correlationPremium;maxDrawdownVal;`OK;"")
 };

.strategy.register[`dispersion;.strategy.dispersion.init;.strategy.dispersion.step;.strategy.dispersion.summary;.strategy.dispersion.defaultConfig];

/ ==================================================================
/ 19. Path adapter: fromFuturesCurve (commodity curve evolution)
/ ==================================================================
/ Evolve a synthetic futures curve over time. evolutionModel:
/   `simple - GBM front + additive contango per tenor.
/   `mrjump - mean-reverting-jump front (via .commodity.mrjump.simulatePaths,
/             pathCount=1) + additive contango per tenor.
/ Returns: frontPath (standard schema; spot=front level), curveSnapshots (long-
/ format table), tenors, evolutionModel, params, jumpCountsAtStep (mrjump only).

.strategy.path.fromFuturesCurve:{[pathCfg]
    required:`tenors`evolutionModel`evolutionParams`steps`stepYears`seed;
    missing:required where not required in key pathCfg;
    if[0<count missing; '"strategy.path.fromFuturesCurve missing keys: ",", " sv string missing];
    tenorsVec:`float$pathCfg`tenors;
    evolModel:pathCfg`evolutionModel;
    evolParams:pathCfg`evolutionParams;
    stepsTotal:pathCfg`steps;
    dt:`float$pathCfg`stepYears;
    seed:pathCfg`seed;
    if[stepsTotal<=1; '"strategy.path.fromFuturesCurve steps must be > 1"];
    if[dt<=0f; '"strategy.path.fromFuturesCurve stepYears must be positive"];
    if[0=count tenorsVec; '"strategy.path.fromFuturesCurve tenors must be non-empty"];
    if[any tenorsVec<=0f; '"strategy.path.fromFuturesCurve tenors must be positive"];
    if[not evolModel in `simple`mrjump; '"strategy.path.fromFuturesCurve evolutionModel must be simple or mrjump"];
    startDate:$[`startDate in key pathCfg; pathCfg`startDate; 2024.01.01];
    stepDates:startDate+til stepsTotal;
    contango:$[`contango in key evolParams; `float$evolParams`contango; 0f];
    jumpCountsAtStep:stepsTotal#0;
    if[evolModel=`simple;
        spot0:`float$evolParams`spot0;
        drift:`float$evolParams`drift;
        vol:`float$evolParams`volatility;
        if[spot0<=0f; '"fromFuturesCurve simple spot0 must be positive"];
        if[vol<0f; '"fromFuturesCurve simple volatility must be non-negative"];
        normals:.montecarlo.__generateNormals[stepsTotal-1;seed];
        logIncrs:((drift-(0.5*vol*vol))*dt)+(vol*sqrt dt)*normals;
        logFronts:(log spot0)+sums logIncrs;
        frontLevels:spot0,exp logFronts];
    if[evolModel=`mrjump;
        x0:`float$evolParams`initialLogPrice;
        mcConfig:`pathCount`timeStepCount`randomSeed`antithetic`confidenceLevel!(1;stepsTotal-1;seed;0b;0.95);
        mrjumpParams:`meanReversionSpeed`longRunLogMean`volatility`jumpIntensity`jumpMean`jumpVolatility#evolParams;
        simResult:.commodity.mrjump.simulatePaths[x0;mrjumpParams;dt*stepsTotal-1;mcConfig];
        pathMat:simResult`pathMatrix;
        jumpMat:simResult`jumpCountMatrix;
        evolvedFront:(pathMat 0);
        frontLevels:(exp x0),evolvedFront;
        jumpCountsAtStep:0i,`int$jumpMat 0];
    rfrConst:$[`riskFreeRate in key evolParams; `float$evolParams`riskFreeRate; 0f];
    frontPath:flip .strategy.path.__schemaCols!(
        til stepsTotal;
        stepDates;
        frontLevels;
        stepsTotal#0Nf;
        stepsTotal#rfrConst;
        stepsTotal#0f;
        stepsTotal#0Nf;
        stepsTotal#`OK);
    / Build long-format curveSnapshots: per (step, tenor) row.
    nTenors:count tenorsVec;
    nRowsLong:stepsTotal*nTenors;
    snapshotStepIdx:nRowsLong#0;
    snapshotDates:nRowsLong#startDate;
    snapshotTenors:nRowsLong#0f;
    snapshotPrices:nRowsLong#0f;
    rowFiller:0;
    stepIdxLoop:0;
    while[stepIdxLoop<stepsTotal;
        tenorIdxLoop:0;
        while[tenorIdxLoop<nTenors;
            snapshotStepIdx[rowFiller]:stepIdxLoop;
            snapshotDates[rowFiller]:stepDates stepIdxLoop;
            snapshotTenors[rowFiller]:tenorsVec tenorIdxLoop;
            snapshotPrices[rowFiller]:(frontLevels stepIdxLoop)+contango*tenorsVec tenorIdxLoop;
            rowFiller+:1;
            tenorIdxLoop+:1];
        stepIdxLoop+:1];
    curveSnapshots:flip `stepIndex`stepDate`tenor`futuresPrice!(snapshotStepIdx;snapshotDates;snapshotTenors;snapshotPrices);
    `frontPath`curveSnapshots`tenors`evolutionModel`evolutionParams`frontLevels`jumpCountsAtStep!(
        frontPath;curveSnapshots;tenorsVec;evolModel;evolParams;frontLevels;jumpCountsAtStep)
 };

/ ==================================================================
/ 20. Concrete strategy: powerSpikeCapture (mrjump-driven)
/ ==================================================================
/ Long calls on the front future (Black-76 priced) when front deviates above
/ mean by deviationThreshold sigmas OR a jump event was seen at this step.
/ Delta-hedge with the front future via __hedgeStep. curveBundle from
/ stratCfg holds the mrjump curve path + jump counts.

.strategy.powerSpikeCapture.__rowEmitCols:`stepIndex`stepDate`spot`frontLevel`callOptionValue`delta`spikeFlag`positionValue`netDelta`hedgePosition`hedgeTrade`txnCost`positionPnl`hedgePnl`financingPnl`stepPnl`cumulativePnl`status`message;

.strategy.powerSpikeCapture.defaultConfig:{[]
    `curveBundle`callStrikePct`callExpiry`callVol`deviationThreshold`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears!(
        ()!();
        0.05;
        0.10;
        0.40;
        2f;
        `interval;1;0.05;0f;0f;1f%252f)
 };

.strategy.powerSpikeCapture.__validateConfig:{[stratCfg]
    requiredKeys:`curveBundle`callStrikePct`callExpiry`callVol`deviationThreshold`rebalanceMode`rebalanceInterval`deltaBand`txnCostRate`financingRate`stepYears;
    missing:requiredKeys where not requiredKeys in key stratCfg;
    if[0<count missing; '"powerSpikeCapture config missing keys: ",", " sv string missing];
    bundle:stratCfg`curveBundle;
    bundleKeys:`frontPath`curveSnapshots`tenors`evolutionModel`frontLevels`jumpCountsAtStep;
    bundleMissing:bundleKeys where not bundleKeys in key bundle;
    if[0<count bundleMissing; '"powerSpikeCapture curveBundle missing keys: ",", " sv string bundleMissing];
    if[(stratCfg`callExpiry)<=0f; '"powerSpikeCapture callExpiry must be positive"];
    if[(stratCfg`callVol)<=0f; '"powerSpikeCapture callVol must be positive"];
    if[(stratCfg`stepYears)<=0f; '"powerSpikeCapture stepYears must be positive"];
 };

.strategy.powerSpikeCapture.__priceCall:{[fwdPrice;strikePrice;expiry;vol;rfr]
    px:.commodity.black76.price[`call;fwdPrice;strikePrice;expiry;vol;rfr];
    greeks:.commodity.black76.greeks[`call;fwdPrice;strikePrice;expiry;vol;rfr];
    `value`delta!(px;greeks`delta)
 };

.strategy.powerSpikeCapture.__isSpike:{[front;mean;sd;jumpCount;threshold]
    deviationSpike:(front-mean)>threshold*sd;
    jumpSpike:jumpCount>0;
    deviationSpike or jumpSpike
 };

.strategy.powerSpikeCapture.init:{[trade;firstStep;model;fdmConfig;stratCfg]
    .strategy.powerSpikeCapture.__validateConfig stratCfg;
    bundle:stratCfg`curveBundle;
    frontLevels:bundle`frontLevels;
    jumpCounts:bundle`jumpCountsAtStep;
    meanFront:avg frontLevels;
    sdFront:dev frontLevels;
    notional:trade`notional;
    spot:firstStep`spot;
    rfr:firstStep`riskFreeRate;
    callStrikePct:stratCfg`callStrikePct;
    callExpiry:stratCfg`callExpiry;
    callVol:stratCfg`callVol;
    strike:spot*1f+callStrikePct;
    spikeFlag:.strategy.powerSpikeCapture.__isSpike[spot;meanFront;sdFront;jumpCounts firstStep`stepIndex;stratCfg`deviationThreshold];
    if[not spikeFlag;
        flatRow:.strategy.powerSpikeCapture.__rowEmitCols!(
            firstStep`stepIndex;firstStep`stepDate;spot;spot;0Nf;0f;0b;
            0f;0f;0f;0f;0f;
            0f;0f;0f;0f;0f;
            `flat;"No entry spike");
        :`notional`callStrikePct`callExpiry`callVol`meanFront`sdFront`deviationThreshold`bundle`spikeCount`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue`prevHeld`hedgeTrade`txnCost`financingPnl`hedgePnl`positionPnl`stepPnl`cumulativePnl`rowEmit!(
            notional;callStrikePct;callExpiry;callVol;meanFront;sdFront;stratCfg`deviationThreshold;bundle;0;0f;0f;0f;0;spot;0f;0b;0f;0f;0f;0f;0f;0f;0f;flatRow)];
    callMark:.strategy.powerSpikeCapture.__priceCall[spot;strike;callExpiry;callVol;rfr];
    callValue:callMark`value;
    callDelta:callMark`delta;
    positionValue:notional*callValue;
    positionDelta:notional*callDelta;
    hedgeInit:.strategy.__hedgeInit `spot`positionDelta`txnCostRate!(spot;positionDelta;stratCfg`txnCostRate);
    initialCash:(neg positionValue)+hedgeInit`cashAdj;
    initialStepPnl:neg hedgeInit`txnCost;
    rowEmit:.strategy.powerSpikeCapture.__rowEmitCols!(
        firstStep`stepIndex;firstStep`stepDate;spot;spot;callValue;callDelta;1b;
        positionValue;positionDelta;hedgeInit`hedgePosition;hedgeInit`hedgeTrade;hedgeInit`txnCost;
        0f;0f;0f;initialStepPnl;initialStepPnl;`OK;"entry spike");
    `notional`callStrikePct`callExpiry`callVol`meanFront`sdFront`deviationThreshold`bundle`spikeCount`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue`prevHeld`hedgeTrade`txnCost`financingPnl`hedgePnl`positionPnl`stepPnl`cumulativePnl`rowEmit!(
        notional;callStrikePct;callExpiry;callVol;meanFront;sdFront;stratCfg`deviationThreshold;bundle;1;initialCash;hedgeInit`hedgePosition;positionDelta;1;spot;positionValue;1b;hedgeInit`hedgeTrade;hedgeInit`txnCost;0f;0f;0f;initialStepPnl;initialStepPnl;rowEmit)
 };

.strategy.powerSpikeCapture.step:{[state;marketStep;trade;model;fdmConfig;stratCfg]
    bundle:state`bundle;
    jumpCounts:bundle`jumpCountsAtStep;
    notional:state`notional;
    spot:marketStep`spot;
    rfr:marketStep`riskFreeRate;
    callStrikePct:state`callStrikePct;
    callExpiry:state`callExpiry;
    callVol:state`callVol;
    stepIndex:marketStep`stepIndex;
    spikeFlag:.strategy.powerSpikeCapture.__isSpike[spot;state`meanFront;state`sdFront;jumpCounts stepIndex;state`deviationThreshold];
    cashPrev:state`cash;
    financingRate:stratCfg`financingRate;
    stepYears:stratCfg`stepYears;
    financingPnl:(financingRate*cashPrev)*stepYears;
    if[(not state`prevHeld) and spikeFlag;
        / Opening on this spike: enter long ATM-ish call, set hedge.
        strike:spot*1f+callStrikePct;
        callMark:.strategy.powerSpikeCapture.__priceCall[spot;strike;callExpiry;callVol;rfr];
        callValue:callMark`value;
        callDelta:callMark`delta;
        newPositionValue:notional*callValue;
        newPositionDelta:notional*callDelta;
        / Open hedge against the new position
        newHedgePos:neg newPositionDelta;
        hedgeTrade:newHedgePos;
        txnCostOpen:((abs hedgeTrade)*spot+notional*callValue)*stratCfg`txnCostRate;
        legCashOpen:neg ((newHedgePos*spot)+notional*callValue);
        newCash:((cashPrev+financingPnl)+legCashOpen)-txnCostOpen;
        / On opening from flat, position-value jump is offset by cash outflow:
        / report positionPnl/hedgePnl as 0 to match PV identity.
        positionPnl:0f;
        hedgePnl:0f;
        stepPnl:financingPnl-txnCostOpen;
        cumulativePnl:(state`cumulativePnl)+stepPnl;
        rowEmit:.strategy.powerSpikeCapture.__rowEmitCols!(
            stepIndex;marketStep`stepDate;spot;spot;callValue;callDelta;1b;
            newPositionValue;newPositionDelta;newHedgePos;hedgeTrade;txnCostOpen;
            positionPnl;hedgePnl;financingPnl;stepPnl;cumulativePnl;`OK;"opened on spike");
        :@[state;`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue`prevHeld`spikeCount`hedgeTrade`txnCost`financingPnl`hedgePnl`positionPnl`stepPnl`cumulativePnl`rowEmit;:;
            (newCash;newHedgePos;newPositionDelta;1+state`numRebalances;spot;newPositionValue;1b;1+state`spikeCount;hedgeTrade;txnCostOpen;financingPnl;hedgePnl;positionPnl;stepPnl;cumulativePnl;rowEmit)]];
    if[state`prevHeld;
        / Already long: mark + hedge step
        strike:(state`prevSpot)*1f+callStrikePct;  / strike fixed at entry
        callMark:.strategy.powerSpikeCapture.__priceCall[spot;strike;callExpiry;callVol;rfr];
        callValue:callMark`value;
        callDelta:callMark`delta;
        newPositionValue:notional*callValue;
        newPositionDelta:notional*callDelta;
        hedgeState:`cash`hedgePosition`hedgedDelta`numRebalances`prevSpot`prevPositionValue#state;
        stepInputs:`spot`positionValue`positionDelta`stepIndex`stepYears`txnCostRate`financingRate`rebalanceMode`rebalanceInterval`deltaBand!(
            spot;newPositionValue;newPositionDelta;stepIndex;stepYears;stratCfg`txnCostRate;financingRate;stratCfg`rebalanceMode;stratCfg`rebalanceInterval;stratCfg`deltaBand);
        newHedge:.strategy.__hedgeStep[hedgeState;stepInputs];
        cumulativePnl:(state`cumulativePnl)+newHedge`stepPnl;
        rowEmit:.strategy.powerSpikeCapture.__rowEmitCols!(
            stepIndex;marketStep`stepDate;spot;spot;callValue;callDelta;spikeFlag;
            newPositionValue;newPositionDelta;newHedge`hedgePosition;newHedge`hedgeTrade;newHedge`txnCost;
            newHedge`positionPnl;newHedge`hedgePnl;newHedge`financingPnl;newHedge`stepPnl;cumulativePnl;`OK;"held");
        :state,newHedge,`cumulativePnl`spikeCount`rowEmit!(cumulativePnl;(state`spikeCount)+`long$spikeFlag;rowEmit)];
    / Else: flat
    newCash:cashPrev+financingPnl;
    stepPnl:financingPnl;
    cumulativePnl:(state`cumulativePnl)+stepPnl;
    rowEmit:.strategy.powerSpikeCapture.__rowEmitCols!(
        stepIndex;marketStep`stepDate;spot;spot;0Nf;0f;0b;
        0f;0f;0f;0f;0f;
        0f;0f;financingPnl;stepPnl;cumulativePnl;`OK;"flat");
    @[state;`cash`prevSpot`prevPositionValue`prevHeld`hedgeTrade`txnCost`financingPnl`hedgePnl`positionPnl`stepPnl`cumulativePnl`rowEmit;:;
        (newCash;spot;0f;0b;0f;0f;financingPnl;0f;0f;stepPnl;cumulativePnl;rowEmit)]
 };

.strategy.powerSpikeCapture.summary:{[resultTable;stratCfg]
    base:`strategyName`steps`spikeCount`spikePremiumCaptured`gateOpen`totalPnl`positionPnlTotal`hedgePnlTotal`financingTotal`txnCostTotal`maxDrawdown`status`errorMessage!(
        `powerSpikeCapture;0;0;0Nf;0b;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;`ERROR;"empty");
    if[(0=count resultTable)|not 98h=type resultTable; :base];
    okRows:resultTable where (resultTable`status)=`OK;
    stepCount:count resultTable;
    if[0=count okRows; :@[base;`steps;:;stepCount]];
    spikeRows:resultTable where resultTable`spikeFlag;
    spikeCountVal:count spikeRows;
    spikePremiumCaptured:$[0<spikeCountVal; sum spikeRows`callOptionValue; 0f];
    totalsDict:first 0!select
        totalPnl:sum stepPnl, positionPnlTotal:sum positionPnl, hedgePnlTotal:sum hedgePnl, financingTotal:sum financingPnl, txnCostTotal:sum txnCost
        from okRows;
    cumPnl:sums okRows`stepPnl;
    maxDrawdownVal:max (maxs cumPnl)-cumPnl;
    gateOpen:0<spikeCountVal;
    totalsDict,`strategyName`steps`spikeCount`spikePremiumCaptured`gateOpen`maxDrawdown`status`errorMessage!(
        `powerSpikeCapture;stepCount;spikeCountVal;spikePremiumCaptured;gateOpen;maxDrawdownVal;`OK;"")
 };

.strategy.register[`powerSpikeCapture;.strategy.powerSpikeCapture.init;.strategy.powerSpikeCapture.step;.strategy.powerSpikeCapture.summary;.strategy.powerSpikeCapture.defaultConfig];

/ ==================================================================
/ 21. Concrete strategy: commodityCalendar (calendar spread on curve)
/ ==================================================================
/ Long 1 futures at nearTenor, short 1 futures at farTenor of the same curve.
/ Each leg's value = spotFutures - entryFutures (mark-to-market). Roll the near
/ leg on expiry (when the step's relative time exceeds nearTenor): close the
/ near at current near-curve price (settled to intrinsic vs entry), open the
/ next tenor as new near. Reuse calendarRoll's leg-table lifecycle + portfolio-
/ value-identity pattern.

.strategy.commodityCalendar.__rowEmitCols:`stepIndex`stepDate`spot`nearFutures`farFutures`nearEntry`farEntry`positionValue`spreadValue`rollEvents`txnCost`positionPnl`rollPnl`financingPnl`stepPnl`cumulativePnl`status`message;

.strategy.commodityCalendar.defaultConfig:{[]
    `curveBundle`nearTenor`farTenor`rollTriggerTenor`txnCostRate`financingRate`stepYears!(
        ()!();
        0.10;
        0.30;
        0.05;
        0f;0f;1f%252f)
 };

.strategy.commodityCalendar.__validateConfig:{[stratCfg]
    requiredKeys:`curveBundle`nearTenor`farTenor`rollTriggerTenor`txnCostRate`financingRate`stepYears;
    missing:requiredKeys where not requiredKeys in key stratCfg;
    if[0<count missing; '"commodityCalendar config missing keys: ",", " sv string missing];
    if[(stratCfg`nearTenor)<=0f; '"commodityCalendar nearTenor must be positive"];
    if[(stratCfg`farTenor)<=stratCfg`nearTenor; '"commodityCalendar farTenor must exceed nearTenor"];
    if[(stratCfg`stepYears)<=0f; '"commodityCalendar stepYears must be positive"];
    bundle:stratCfg`curveBundle;
    bundleKeys:`frontPath`curveSnapshots`tenors`evolutionParams`frontLevels;
    bundleMissing:bundleKeys where not bundleKeys in key bundle;
    if[0<count bundleMissing; '"commodityCalendar curveBundle missing keys: ",", " sv string bundleMissing];
 };

.strategy.commodityCalendar.__priceAtTenor:{[bundle;stepIndex;targetTenor]
    snapshots:bundle`curveSnapshots;
    stepRows:snapshots where (snapshots`stepIndex)=stepIndex;
    if[0=count stepRows; :0Nf];
    tenorsAvail:stepRows`tenor;
    pricesAvail:stepRows`futuresPrice;
    diffs:abs tenorsAvail-targetTenor;
    bestIdx:diffs?min diffs;
    pricesAvail bestIdx
 };

.strategy.commodityCalendar.init:{[trade;firstStep;model;fdmConfig;stratCfg]
    .strategy.commodityCalendar.__validateConfig stratCfg;
    bundle:stratCfg`curveBundle;
    notional:trade`notional;
    nearTenor:stratCfg`nearTenor;
    farTenor:stratCfg`farTenor;
    stepIndex0:firstStep`stepIndex;
    nearEntry:.strategy.commodityCalendar.__priceAtTenor[bundle;stepIndex0;nearTenor];
    farEntry:.strategy.commodityCalendar.__priceAtTenor[bundle;stepIndex0;farTenor];
    spot:firstStep`spot;
    legBookInit:`name`tenor`entryPrice`units`active!(
        `near`far;nearTenor,farTenor;nearEntry,farEntry;notional,neg notional;1b,1b);
    legBookInit:flip legBookInit;
    / At entry: each leg's mark = current - entry = 0
    positionValue:0f;
    spreadValueEntry:farEntry-nearEntry;
    legEntryTxnCost:((abs notional*nearEntry)+(abs notional*farEntry))*stratCfg`txnCostRate;
    initialStepPnl:neg legEntryTxnCost;
    initialCash:neg legEntryTxnCost;  / futures positions have zero cash flow at entry; only txn costs
    rowEmit:.strategy.commodityCalendar.__rowEmitCols!(
        stepIndex0;firstStep`stepDate;spot;nearEntry;farEntry;nearEntry;farEntry;
        positionValue;spreadValueEntry;0;legEntryTxnCost;
        0f;0f;0f;initialStepPnl;initialStepPnl;`OK;"");
    `notional`nearTenor`farTenor`rollTriggerTenor`bundle`legBook`rollEventCount`cash`prevSpot`prevPositionValue`hedgePosition`txnCost`financingPnl`positionPnl`rollPnl`stepPnl`cumulativePnl`rowEmit!(
        notional;nearTenor;farTenor;stratCfg`rollTriggerTenor;bundle;legBookInit;0;initialCash;spot;positionValue;0f;legEntryTxnCost;0f;0f;0f;initialStepPnl;initialStepPnl;rowEmit)
 };

.strategy.commodityCalendar.step:{[state;marketStep;trade;model;fdmConfig;stratCfg]
    bundle:state`bundle;
    legBook:state`legBook;
    notional:state`notional;
    nearTenor:state`nearTenor;
    farTenor:state`farTenor;
    stepIndex:marketStep`stepIndex;
    spot:marketStep`spot;
    cashPrev:state`cash;
    financingRate:stratCfg`financingRate;
    stepYears:stratCfg`stepYears;
    txnCostRate:stratCfg`txnCostRate;
    / Step-time-elapsed for the front (near) tenor concept
    elapsedYears:stepIndex*stepYears;
    / Roll near leg if it has expired (i.e. elapsedYears >= nearTenor when we started at near)
    nearLegRow:first 0!legBook where (legBook`name)=`near;
    farLegRow:first 0!legBook where (legBook`name)=`far;
    nearShouldRoll:elapsedYears>=nearLegRow`tenor;
    rollPnl:0f;
    rollTxnCost:0f;
    rollEventsCount:0;
    nearEntryEff:nearLegRow`entryPrice;
    if[nearShouldRoll;
        / Settle near at near-curve price (this step's price for that tenor) -> close near leg
        nearSettlePrice:.strategy.commodityCalendar.__priceAtTenor[bundle;stepIndex;nearTenor];
        nearLegMark:(nearLegRow`units)*nearSettlePrice-nearEntryEff;
        rollPnl+:nearLegMark;
        / Open a new near at new tenor = (state`rollTriggerTenor will be used for next length)
        nextNearTenor:state`rollTriggerTenor;
        newNearEntry:.strategy.commodityCalendar.__priceAtTenor[bundle;stepIndex;nextNearTenor];
        rollTxnCost+:(abs nearLegRow`units)*nearSettlePrice*txnCostRate;
        rollTxnCost+:(abs notional)*newNearEntry*txnCostRate;
        legBook:update entryPrice:newNearEntry, tenor:nextNearTenor from legBook where name=`near;
        nearEntryEff:newNearEntry;
        rollEventsCount:1];
    / Re-fetch nearLegRow after possible roll
    nearLegRow:first 0!legBook where (legBook`name)=`near;
    farLegRow:first 0!legBook where (legBook`name)=`far;
    nearMarkPrice:.strategy.commodityCalendar.__priceAtTenor[bundle;stepIndex;nearLegRow`tenor];
    farMarkPrice:.strategy.commodityCalendar.__priceAtTenor[bundle;stepIndex;farLegRow`tenor];
    nearLegValue:(nearLegRow`units)*nearMarkPrice-nearLegRow`entryPrice;
    farLegValue:(farLegRow`units)*farMarkPrice-farLegRow`entryPrice;
    newPositionValue:nearLegValue+farLegValue;
    spreadValue:farMarkPrice-nearMarkPrice;
    positionPnl:(newPositionValue-state`prevPositionValue)-rollPnl;
    financingPnl:(financingRate*cashPrev)*stepYears;
    / Cash flow at roll = nearLegMark (received) - newNearEntry * 0 (futures are zero-cost at open) - rollTxnCost
    cashFromRoll:rollPnl-rollTxnCost;
    newCash:(cashPrev+financingPnl)+cashFromRoll;
    txnCostStep:rollTxnCost;
    stepPnl:((positionPnl+rollPnl)+financingPnl)-txnCostStep;
    cumulativePnl:(state`cumulativePnl)+stepPnl;
    rowEmit:.strategy.commodityCalendar.__rowEmitCols!(
        stepIndex;marketStep`stepDate;spot;nearMarkPrice;farMarkPrice;nearLegRow`entryPrice;farLegRow`entryPrice;
        newPositionValue;spreadValue;rollEventsCount;txnCostStep;
        positionPnl;rollPnl;financingPnl;stepPnl;cumulativePnl;`OK;"");
    @[state;`legBook`rollEventCount`cash`prevSpot`prevPositionValue`hedgePosition`txnCost`financingPnl`positionPnl`rollPnl`stepPnl`cumulativePnl`rowEmit;:;
        (legBook;(state`rollEventCount)+rollEventsCount;newCash;spot;newPositionValue;0f;txnCostStep;financingPnl;positionPnl;rollPnl;stepPnl;cumulativePnl;rowEmit)]
 };

.strategy.commodityCalendar.summary:{[resultTable;stratCfg]
    base:`strategyName`steps`nearTenor`farTenor`totalRolls`spreadPnlTotal`rollPnlTotal`totalPnl`positionPnlTotal`financingTotal`txnCostTotal`maxDrawdown`status`errorMessage!(
        `commodityCalendar;0;stratCfg`nearTenor;stratCfg`farTenor;0;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;`ERROR;"empty");
    if[(0=count resultTable)|not 98h=type resultTable; :base];
    okRows:resultTable where (resultTable`status)=`OK;
    stepCount:count resultTable;
    if[0=count okRows; :@[base;`steps;:;stepCount]];
    totalsDict:first 0!select
        totalPnl:sum stepPnl, positionPnlTotal:sum positionPnl, rollPnlTotal:sum rollPnl, financingTotal:sum financingPnl, txnCostTotal:sum txnCost
        from okRows;
    totalRollsVal:sum okRows`rollEvents;
    spreadPnlTotal:totalsDict`positionPnlTotal;
    cumPnl:sums okRows`stepPnl;
    maxDrawdownVal:max (maxs cumPnl)-cumPnl;
    totalsDict,`strategyName`steps`nearTenor`farTenor`totalRolls`spreadPnlTotal`maxDrawdown`status`errorMessage!(
        `commodityCalendar;stepCount;stratCfg`nearTenor;stratCfg`farTenor;totalRollsVal;spreadPnlTotal;maxDrawdownVal;`OK;"")
 };

.strategy.register[`commodityCalendar;.strategy.commodityCalendar.init;.strategy.commodityCalendar.step;.strategy.commodityCalendar.summary;.strategy.commodityCalendar.defaultConfig];
