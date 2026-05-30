/ commodityStrategies.q - real-data commodity strategy backtests (v0.53)
/ ----------------------------------------------------------------------------
/ Tests whether the calibration / Kalman outputs are tradeable SIGNALS by
/ backtesting commodity strategies on the real WTI curve history and ranking
/ them. Strategies self-register on the existing generic engine (lib/strategy.q)
/ and read PRECOMPUTED, CAUSAL signals from a signal-augmented path; no model is
/ ever re-fit inside the backtest loop.
/ ----------------------------------------------------------------------------
/ ASSUMPTIONS (documented):
/   ROLL CONVENTION: the tradeable continuous series is ROLL-ADJUSTED RETURNS -
/     hold the front contract, roll to the next a configurable number of calendar
/     days before expiry; the daily return is always the held contract's own
/     return (no level-stitch jump). The roll yield is captured in the returns.
/   CAUSALITY / NO LOOK-AHEAD: every signal at date t uses only data up to t. The
/     Kalman params are estimated ONCE on a TRAIN split (fixed trainEndDate so the
/     split does not move when future data is appended) and the causal forward
/     filter produces chi/xi; per-date convenience yield is an independent
/     single-curve fit; momentum/carry use only trailing/contemporaneous data.
/   VOL TARGETING: positions are scaled by targetVol / (annualised train-window
/     front-return vol) so all strategies trade at a common vol - the ranking
/     compares edge, not leverage. Positions are LAGGED one step: stepPnl_t =
/     position_{t-1} * frontReturn_t - turnoverCost_t.
/   COSTS: turnover cost = |dPosition| * txnCostRate each step.
/   FUTURES MTM: P&L is realised to cash each step; the futures position carries
/     zero mark, so PV = cash and deltaPV == stepPnl (the accounting identity).
/ ----------------------------------------------------------------------------
/ Reuses (never modifies): .parser.crude.*, .commodity.curveCal.convenienceYield-
/   Series, .commodity.kalman.{panelFromCurveHistory,estimate,filter},
/   .strategy.{run,register,__portfolioValue}, commodityCalendar.

/ ============================================================================
/ A0. Signal-augmented real-curve path builder
/ ============================================================================

.strategy.path.__commodityDefaultSigCfg:{[]
    `rollDaysBeforeExpiry`trainFraction`trainEndDate`momentumLookback`txnCostRate`targetVol`riskFreeRate`storageCostRate`annualizationDays`kalmanEstCfg`carryMargin!(
        5;
        0.6;
        0Nd;
        20;
        0.0005;
        0.15;
        0.02;
        0.01;
        252f;
        `gridSteps`refineRounds`nSweeps!(7;3;3);
        0.0)
 };

/ Per-date contractYM!price dictionary for fast held-contract lookups.
.strategy.path.__priceMapForDate:{[curveHist;asofVal]
    sub:select contractYM,price from curveHist where asofDate=asofVal;
    (sub`contractYM)!sub`price
 };

/ Roll decision: given the contract held into day dayIdx, return the contract to
/ hold into the NEXT day (roll to the next tenor when within rollDays of expiry).
.strategy.path.__rollNext:{[perDateList;datesVec;rollDays;heldYM;dayIdx]
    dayTbl:perDateList dayIdx;
    heldRow:select from dayTbl where contractYM=heldYM;
    if[0=count heldRow; :heldYM];
    daysToExpiry:(first heldRow`expiry)-datesVec dayIdx;
    if[daysToExpiry>rollDays; :heldYM];
    laterContracts:select from dayTbl where tenor>first heldRow`tenor;
    if[0=count laterContracts; :heldYM];
    first (`tenor xasc laterContracts)`contractYM
 };

.strategy.path.commoditySignals:{[curveHistory;sigCfg]
    if[not 98h=type curveHistory; '"commoditySignals: curveHistory must be a table"];
    if[not all `asofDate`tenor`price`contractYM`expiry in cols curveHistory;
        '"commoditySignals: curveHistory needs asofDate, tenor, price, contractYM, expiry"];
    cfg:.strategy.path.__commodityDefaultSigCfg[];
    if[count sigCfg; cfg:cfg,sigCfg];
    / Exclude non-positive prices (log domain) and sort.
    curveHist:`asofDate`tenor xasc select from curveHistory where price>0f;
    dates:asc distinct curveHist`asofDate;
    nDates:count dates;
    if[nDates<3; '"commoditySignals: need >=3 dates"];
    / train/test split: fixed trainEndDate (causally stable) or derived from fraction.
    trainEndDate:$[not null cfg`trainEndDate; cfg`trainEndDate; dates floor (cfg`trainFraction)*nDates-1];
    isTrain:dates<=trainEndDate;
    / Per-date tenor-sorted contract tables.
    perDateList:{[ch;d] `tenor xasc select contractYM,price,tenor,expiry from ch where asofDate=d}[curveHist;] each dates;
    priceMap:dates!.strategy.path.__priceMapForDate[curveHist;] each dates;
    / --- roll-adjusted front returns (causal walk, runs once) ---
    front0:first (perDateList 0)`contractYM;
    heldTail:.strategy.path.__rollNext[perDateList;dates;cfg`rollDaysBeforeExpiry;;]\[front0;til nDates-1];
    heldSeq:front0,heldTail;
    frontReturn:0n,{[priceMap;dates;heldSeq;t]
        ((priceMap[dates t]heldSeq t)%priceMap[dates t-1]heldSeq t)-1f}[priceMap;dates;heldSeq;] each 1+til nDates-1;
    frontPrice:{[priceMap;dates;heldSeq;t] priceMap[dates t]heldSeq t}[priceMap;dates;heldSeq;] each til nDates;
    daysToExpiry:`float${[perDateList;dates;heldSeq;t]
        hr:select from perDateList[t] where contractYM=heldSeq t;
        $[0=count hr;0Ni;(first hr`expiry)-dates t]}[perDateList;dates;heldSeq;] each til nDates;
    / --- model-free curve carry: annualised near-vs-deferred backwardation ---
    curveSlopeCarry:{[perDateList;t]
        pd:perDateList t;
        if[2>count pd;:0n];
        nearP:pd[`price]0; farP:pd[`price]1; nearTau:pd[`tenor]0; farTau:pd[`tenor]1;
        ((log nearP)-log farP)%farTau-nearTau}[perDateList;] each til nDates;
    / --- momentum: trailing-N mean daily return (causal) ---
    momentum:(cfg`momentumLookback) mavg 0f^frontReturn;
    / --- Kalman: estimate on TRAIN, filter forward on FULL panel (causal) ---
    panel:.commodity.kalman.panelFromCurveHistory curveHist;
    trainPanel:select from panel where obsDate<=trainEndDate;
    kalEst:.commodity.kalman.estimate[trainPanel;cfg`kalmanEstCfg];
    kalParams:kalEst`estimatedParams;
    kalFwd:.commodity.kalman.filter[panel;kalParams];
    / align filter output (its dates) to our date axis
    filterByDate:(kalFwd`dates)!flip (kalFwd`chi;kalFwd`xi);
    chiXi:filterByDate dates;
    chi:chiXi[;0]; xi:chiXi[;1];
    kappaEst:kalParams`kappa; sigChiEst:kalParams`sigChi;
    stationaryStd:sigChiEst%sqrt 2f*kappaEst;
    chiZ:$[stationaryStd>0f; chi%stationaryStd; nDates#0f];
    / --- convenience yield per date (single-curve fit, uses estimated kappa) ---
    cyCalCfg:`kappa`shortVolatility`longVolatility`correlation`riskFreeRate!(
        kappaEst;sigChiEst;kalParams`sigXi;kalParams`correlation;cfg`riskFreeRate);
    cyOut:.commodity.curveCal.convenienceYieldSeries[curveHist;cyCalCfg];
    cySeries:cyOut`series;
    cyByDate:(cySeries`asofDate)!cySeries`netConvenienceYield;
    convenienceYield:cyByDate[dates];
    / --- vol-target scale from TRAIN front-return vol (annualised) ---
    trainReturns:(0f^frontReturn) where isTrain;
    dailyVol:dev 1_trainReturns;
    annVol:dailyVol*sqrt cfg`annualizationDays;
    volTargetScale:$[annVol>0f; (cfg`targetVol)%annVol; 1f];
    / --- assemble standard-schema path + signal columns ---
    pathTbl:flip `stepIndex`stepDate`spot`volatility`riskFreeRate`dividendYield`marketPrice`status`frontReturn`frontPrice`frontContractYM`daysToExpiry`curveSlopeCarry`momentum`chi`xi`chiZ`convenienceYield`isTrain`volTargetScale!(
        til nDates;
        dates;
        frontPrice;
        nDates#0Nf;
        nDates#`float$cfg`riskFreeRate;
        nDates#0f;
        frontPrice;
        nDates#`OK;
        0f^frontReturn;
        frontPrice;
        heldSeq;
        daysToExpiry;
        0f^curveSlopeCarry;
        momentum;
        chi;
        xi;
        chiZ;
        convenienceYield;
        isTrain;
        nDates#volTargetScale);
    `path`kalmanParams`trainEndDate`nTrain`nTest`volTargetScale`stationaryStd`heldContracts!(
        pathTbl;kalParams;trainEndDate;sum isTrain;sum not isTrain;volTargetScale;stationaryStd;distinct heldSeq)
 };

/ ============================================================================
/ Shared return-backtest core (vol-targeted position * return P&L, futures MTM)
/ ============================================================================
/ Strategies map a signal to a raw target in {-1,0,+1} (or continuous); the core
/ scales by volTargetScale*notional, lags the position one step, charges turnover
/ cost, and realises P&L to cash. PV=cash so deltaPV==stepPnl. Strategies differ
/ only in the raw-target rule; the accounting is shared and identical.

.strategy.commodityBT.__rowEmitCols:`stepIndex`stepDate`frontPrice`signal`rawTarget`position`frontReturn`positionPnl`turnoverCost`stepPnl`cumulativePnl`isTrain`status`message;

.strategy.commodityBT.coreInit:{[trade;firstStep;stratCfg;rawTarget0;signalVal0]
    notional:trade`notional;
    scale:firstStep`volTargetScale;
    txnRate:stratCfg`txnCostRate;
    position0:rawTarget0*scale*notional;
    entryCost:(abs position0)*txnRate;
    stepPnl0:neg entryCost;
    rowEmit:.strategy.commodityBT.__rowEmitCols!(
        firstStep`stepIndex;firstStep`stepDate;firstStep`frontPrice;signalVal0;rawTarget0;position0;
        firstStep`frontReturn;0f;entryCost;stepPnl0;stepPnl0;firstStep`isTrain;`OK;"");
    `cash`prevPosition`prevRawTarget`cumulativePnl`notional`rowEmit!(
        neg entryCost;position0;rawTarget0;stepPnl0;notional;rowEmit)
 };

.strategy.commodityBT.coreStep:{[state;marketStep;stratCfg;rawTarget;signalVal]
    notional:state`notional;
    scale:marketStep`volTargetScale;
    txnRate:stratCfg`txnCostRate;
    prevPos:state`prevPosition;
    ret:marketStep`frontReturn;
    desired:rawTarget*scale*notional;
    positionPnl:prevPos*ret;
    turnover:(abs desired-prevPos)*txnRate;
    stepPnl:positionPnl-turnover;
    newCash:(state`cash)+stepPnl;
    cumulativePnl:(state`cumulativePnl)+stepPnl;
    rowEmit:.strategy.commodityBT.__rowEmitCols!(
        marketStep`stepIndex;marketStep`stepDate;marketStep`frontPrice;signalVal;rawTarget;desired;
        ret;positionPnl;turnover;stepPnl;cumulativePnl;marketStep`isTrain;`OK;"");
    @[state;`cash`prevPosition`prevRawTarget`cumulativePnl`rowEmit;:;(newCash;desired;rawTarget;cumulativePnl;rowEmit)]
 };

/ Time-series performance summary (split into in-sample / out-of-sample by isTrain).
.strategy.commodityBT.coreSummary:{[resultTable;stratCfg;strategyName]
    annDays:$[`annualizationDays in key stratCfg; stratCfg`annualizationDays; 252f];
    notional:$[`notional in key stratCfg; stratCfg`notional; 1f];
    base:`strategyName`steps`totalPnl`testPnl`testAnnualReturn`testAnnualVol`testSharpe`testMaxDrawdown`testHitRate`trainSharpe`turnoverCostTotal`status`errorMessage!(
        strategyName;0;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;`ERROR;"empty");
    if[(0=count resultTable)|not 98h=type resultTable; :base];
    okRows:resultTable where (resultTable`status)=`OK;
    if[0=count okRows; :@[base;`steps;:;count resultTable]];
    perf:.strategy.commodityBT.__perf[;annDays;notional];
    testRows:okRows where not okRows`isTrain;
    trainRows:okRows where okRows`isTrain;
    testPerf:perf testRows`stepPnl;
    trainPerf:perf trainRows`stepPnl;
    `strategyName`steps`totalPnl`testPnl`testAnnualReturn`testAnnualVol`testSharpe`testMaxDrawdown`testHitRate`trainSharpe`turnoverCostTotal`status`errorMessage!(
        strategyName;count okRows;sum okRows`stepPnl;testPerf`totalPnl;testPerf`annualReturn;testPerf`annualVol;
        testPerf`sharpe;testPerf`maxDrawdown;testPerf`hitRate;trainPerf`sharpe;sum okRows`turnoverCost;`OK;"")
 };

/ Time-series stats of a daily P&L vector (returns expressed per unit notional).
.strategy.commodityBT.__perf:{[stepPnlVec;annDays;notional]
    n:count stepPnlVec;
    if[0=n; :`totalPnl`annualReturn`annualVol`sharpe`maxDrawdown`hitRate!(0f;0Nf;0Nf;0Nf;0Nf;0Nf)];
    dailyRet:stepPnlVec%notional;
    meanRet:avg dailyRet;
    volRet:dev dailyRet;
    annualReturn:meanRet*annDays;
    annualVol:volRet*sqrt annDays;
    sharpe:$[volRet>0f; annualReturn%annualVol; 0Nf];
    cumPnl:sums stepPnlVec;
    maxDrawdown:max (maxs cumPnl)-cumPnl;
    hitRate:`float$(sum stepPnlVec>0f)%n;
    `totalPnl`annualReturn`annualVol`sharpe`maxDrawdown`hitRate!(sum stepPnlVec;annualReturn;annualVol;sharpe;maxDrawdown;hitRate)
 };

/ ============================================================================
/ A1. convenienceYieldCarry - long backwardation, flat/short contango
/ ============================================================================
.strategy.convenienceYieldCarry.defaultConfig:{[]
    `signalSource`carryMargin`riskFreeRate`allowShort`txnCostRate`annualizationDays`notional!(
        `convenienceYield;0.0;0.02;1b;0.0005;252f;1f)
 };

.strategy.convenienceYieldCarry.__rawTarget:{[signalVal;stratCfg]
    rate:stratCfg`riskFreeRate; margin:stratCfg`carryMargin;
    $[signalVal>rate+margin; 1f;
      (signalVal<rate-margin) and stratCfg`allowShort; -1f;
      signalVal<rate-margin; 0f;
      0f]
 };

.strategy.convenienceYieldCarry.__signal:{[stepRow;stratCfg] stepRow stratCfg`signalSource};

.strategy.convenienceYieldCarry.init:{[trade;firstStep;model;fdmConfig;stratCfg]
    sig:.strategy.convenienceYieldCarry.__signal[firstStep;stratCfg];
    .strategy.commodityBT.coreInit[trade;firstStep;stratCfg;.strategy.convenienceYieldCarry.__rawTarget[sig;stratCfg];sig]
 };

.strategy.convenienceYieldCarry.step:{[state;marketStep;trade;model;fdmConfig;stratCfg]
    sig:.strategy.convenienceYieldCarry.__signal[marketStep;stratCfg];
    .strategy.commodityBT.coreStep[state;marketStep;stratCfg;.strategy.convenienceYieldCarry.__rawTarget[sig;stratCfg];sig]
 };

.strategy.convenienceYieldCarry.summary:{[resultTable;stratCfg]
    .strategy.commodityBT.coreSummary[resultTable;stratCfg;`convenienceYieldCarry]
 };

.strategy.register[`convenienceYieldCarry;.strategy.convenienceYieldCarry.init;.strategy.convenienceYieldCarry.step;.strategy.convenienceYieldCarry.summary;.strategy.convenienceYieldCarry.defaultConfig];

/ ============================================================================
/ A2. chiReversion - mean-revert the Kalman short-term factor (z = chi/stationaryStd)
/ ============================================================================
.strategy.chiReversion.defaultConfig:{[]
    `entryZ`exitZ`txnCostRate`annualizationDays`notional!(1.0;0.3;0.0005;252f;1f)
 };

/ Hysteresis: enter long when z below -entryZ, short above entryZ, exit inside
/ exitZ, otherwise hold the previous target. (Mean reversion: long a cheap chi.)
.strategy.chiReversion.__rawTarget:{[z;prevRaw;stratCfg]
    entryZ:stratCfg`entryZ; exitZ:stratCfg`exitZ;
    $[z< neg entryZ; 1f;
      z>entryZ; -1f;
      (abs z)<exitZ; 0f;
      prevRaw]
 };

.strategy.chiReversion.init:{[trade;firstStep;model;fdmConfig;stratCfg]
    z:firstStep`chiZ;
    raw:.strategy.chiReversion.__rawTarget[z;0f;stratCfg];
    .strategy.commodityBT.coreInit[trade;firstStep;stratCfg;raw;z]
 };

.strategy.chiReversion.step:{[state;marketStep;trade;model;fdmConfig;stratCfg]
    z:marketStep`chiZ;
    raw:.strategy.chiReversion.__rawTarget[z;state`prevRawTarget;stratCfg];
    .strategy.commodityBT.coreStep[state;marketStep;stratCfg;raw;z]
 };

.strategy.chiReversion.summary:{[resultTable;stratCfg]
    .strategy.commodityBT.coreSummary[resultTable;stratCfg;`chiReversion]
 };

.strategy.register[`chiReversion;.strategy.chiReversion.init;.strategy.chiReversion.step;.strategy.chiReversion.summary;.strategy.chiReversion.defaultConfig];

/ ============================================================================
/ A3. realCurveCalendarRoll - drive the EXISTING commodityCalendar over the real
/ curve history. Minimal new code: a curve-bundle builder from the parser output
/ + a driver that reports realized roll yield split by regime. Accounting is
/ inherited from commodityCalendar (unchanged, already independently tested).
/ ============================================================================

/ Build a commodityCalendar-compatible curve bundle from a parser curveHistory.
.strategy.path.curveBundleFromHistory:{[curveHistory]
    if[not 98h=type curveHistory; '"curveBundleFromHistory: table expected"];
    ch:`asofDate`tenor xasc select from curveHistory where price>0f;
    dates:asc distinct ch`asofDate;
    n:count dates;
    if[n<2; '"curveBundleFromHistory: need >=2 dates"];
    dateIdx:dates!til n;
    curveSnapshots:select stepIndex:dateIdx asofDate, stepDate:asofDate, tenor:`float$tenor, futuresPrice:`float$price from ch;
    / ch is sorted by asofDate then tenor, so the first price per date is the front.
    frontLevels:`float$exec frontP from 0!select frontP:first price by asofDate from ch;
    frontPath:flip .strategy.path.__schemaCols!(
        til n; dates; frontLevels; n#0Nf; n#0.02; n#0f; n#0Nf; n#`OK);
    evolParams:`spot0`drift`volatility`contango`riskFreeRate!(first frontLevels;0f;0.2;0f;0.02);
    `frontPath`curveSnapshots`tenors`evolutionModel`evolutionParams`frontLevels`jumpCountsAtStep!(
        frontPath;curveSnapshots;asc distinct `float$ch`tenor;`real;evolParams;frontLevels;n#0)
 };

.strategy.realCurveCalendarRoll:{[curveHistory;calCfg]
    bundle:.strategy.path.curveBundleFromHistory curveHistory;
    cfg:.strategy.defaultConfig `commodityCalendar;
    if[count calCfg; cfg:cfg,calCfg];
    cfg:@[cfg;`curveBundle;:;bundle];
    trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
        `RCC;`WTI;`equityOption;`european;`call;60f;0.5;1f);
    bsModel:.model.createBlackScholesModel[];
    fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
        `crankNicolson;40;50;0f;300f;`linear;1b;1b);
    runB:.strategy.runAndSummarize[`commodityCalendar;trade;bundle`frontPath;bsModel;fdmCfg;cfg];
    res:runB`result;
    okRes:res where (res`status)=`OK;
    regimeCol:?[okRes[`spreadValue]<0f;`backwardation;`contango];
    rollByRegime:0!select rollPnlTotal:sum rollPnl, positionPnlTotal:sum positionPnl, steps:count i by regime:regimeCol from okRes;
    `result`summary`rollByRegime`bundle!(res;runB`summary;rollByRegime;bundle)
 };

/ ============================================================================
/ Part C. Scoring & ranking on a single real history (time-series metrics).
/ ----------------------------------------------------------------------------
/ The v0.45 ensemble dashboard computes CROSS-PATH stats (mean/std of totalPnl
/ over many synthetic paths) - on one real history that std is 0 and Sharpe is
/ undefined. So single-history ranking needs TIME-SERIES metrics (annualised
/ Sharpe/drawdown/hit-rate from each strategy's daily out-of-sample stepPnl),
/ which the coreSummary already produces; this scorer collects them into a
/ ranked table and a cross-strategy correlation of out-of-sample daily P&L.
/ ============================================================================
.strategy.commodityBT.runSuite:{[strategyNames;trade;path;model;fdmConfig;cfgByName]
    runOne:{[nm;tradeL;pathL;modelL;fdmL;cfgByNameL]
        cfg:$[nm in key cfgByNameL; cfgByNameL nm; .strategy.defaultConfig nm];
        runB:.strategy.runAndSummarize[nm;tradeL;pathL;modelL;fdmL;cfg];
        res:runB`result;
        testPnl:exec stepPnl from res where status=`OK, not isTrain;
        `name`summary`testPnl!(nm;runB`summary;testPnl)
        }[;trade;path;model;fdmConfig;cfgByName];
    runs:runOne each strategyNames;
    perfTbl:.strategy.__rowDictsToTable runs[;`summary];
    ranked:`testSharpe xdesc perfTbl;
    names:runs[;`name];
    pnls:runs[;`testPnl];
    nS:count names;
    corMat:{[i;pnlsL] {[vi;vj] $[(0=dev vi)|0=dev vj; $[vi~vj;1f;0n]; vi cor vj]}[pnlsL i;] each pnlsL}[;pnls] each til nS;
    `ranked`correlationNames`correlationMatrix`results!(ranked;names;corMat;runs)
 };
