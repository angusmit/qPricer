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
