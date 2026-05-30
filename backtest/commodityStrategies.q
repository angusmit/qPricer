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
    .cfg.strategy.commoditySignals
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

/ ---- gas deseasonalization (v0.55) -------------------------------------------
/ Seasonal products (NG) carry a deterministic delivery-month premium that makes
/ the raw curve slope/carry mostly mechanical. We FIT 12 monthly log-premium
/ factors on the TRAIN window only (causal), then divide each contract's price by
/ exp(factor[deliveryMonth]) before extracting the MODEL signals (curve carry,
/ convenience yield, Kalman chi/xi). The tradeable front-return series stays RAW.
/ deseasonalize defaults off (crude) so WTI signal output is byte-identical.
.strategy.path.__monthIdxFromYM:{[contractYM] (contractYM mod 100)-1};

.strategy.path.__fitSeasonalTrain:{[curveHist;trainEndDate]
    trainRows:select from curveHist where asofDate<=trainEndDate, price>0f;
    trainRows:update demeaned:(log price)-avg log price by asofDate from trainRows;
    fracs:(`float$.strategy.path.__monthIdxFromYM trainRows`contractYM)%12f;
    0f^.commodity.seasonality.fitMonthlyFactors[fracs;trainRows`demeaned]
 };

.strategy.path.__deseasonalizeCurve:{[curveHist;monthFactors]
    update price:price%exp monthFactors .strategy.path.__monthIdxFromYM contractYM from curveHist
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
    / Deseasonalised SIGNAL curve (train-only seasonal fit); the RAW curve still
    / drives the tradeable front returns. deseasonalize off -> signal == raw curve.
    deseasonalize:cfg`deseasonalize;
    monthFactors:$[deseasonalize; .strategy.path.__fitSeasonalTrain[curveHist;trainEndDate]; 12#0f];
    signalCurveHist:$[deseasonalize; .strategy.path.__deseasonalizeCurve[curveHist;monthFactors]; curveHist];
    signalPerDateList:$[deseasonalize;
        {[ch;d] `tenor xasc select contractYM,price,tenor,expiry from ch where asofDate=d}[signalCurveHist;] each dates;
        perDateList];
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
    / (computed on the SIGNAL curve - deseasonalised for gas, raw for crude).
    curveSlopeCarry:{[signalPerDateList;t]
        pd:signalPerDateList t;
        if[2>count pd;:0n];
        nearP:pd[`price]0; farP:pd[`price]1; nearTau:pd[`tenor]0; farTau:pd[`tenor]1;
        ((log nearP)-log farP)%farTau-nearTau}[signalPerDateList;] each til nDates;
    / --- momentum: trailing-N mean daily return (causal) ---
    momentum:(cfg`momentumLookback) mavg 0f^frontReturn;
    / --- Kalman: estimate on TRAIN, filter forward on FULL panel (causal) ---
    panel:.commodity.kalman.panelFromCurveHistory signalCurveHist;
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
    cyOut:.commodity.curveCal.convenienceYieldSeries[signalCurveHist;cyCalCfg];
    cySeries:cyOut`series;
    cyByDate:(cySeries`asofDate)!cySeries`netConvenienceYield;
    convenienceYield:cyByDate[dates];
    / --- xi (equilibrium-level) momentum: trailing-N change in the slow factor ---
    xiMomentum:0f^xi-(cfg`momentumLookback) xprev xi;
    / --- deferred (next) contract roll-adjusted return + near-far calendar spread ---
    allContracts:asc distinct curveHist`contractYM;
    nextContract:allContracts!(1_allContracts),0N;
    deferredSeq:nextContract heldSeq;
    deferredReturn:0n,{[priceMap;dates;deferredSeq;t]
        dy:deferredSeq t;
        pNow:priceMap[dates t]dy; pPrev:priceMap[dates t-1]dy;
        $[(null dy)|(null pNow)|null pPrev; 0n; (pNow%pPrev)-1f]}[priceMap;dates;deferredSeq;] each 1+til nDates-1;
    nearFarSpreadReturn:(0f^frontReturn)-0f^deferredReturn;
    / --- vol-target scales from TRAIN-window vol of each traded series (annualised) ---
    annScale:sqrt cfg`annualizationDays;
    volScaleOf:{[retVec;isTrain;targetVol;annScale]
        v:(dev 1_retVec where isTrain)*annScale; $[v>0f; targetVol%v; 1f]}[;isTrain;cfg`targetVol;annScale];
    volTargetScale:volScaleOf 0f^frontReturn;
    volTargetScaleNearFar:volScaleOf 0f^nearFarSpreadReturn;
    / --- curve relative value: per-date schwartz2 residuals -> cheap/rich legs ---
    / error = modelPrice - marketPrice; >0 cheap (long), <0 rich (short). kappa fixed.
    rvCalCfg:cyCalCfg,`meanReversionSpeedRange`gridSteps`refineRounds!((kappaEst;kappaEst);1;1);
    rvIdent:{[curveHistL;rvCalCfgL;d]
        sub:`tenor xasc select tenor,price,contractYM from curveHistL where asofDate=d;
        if[3>count sub; :`longYM`shortYM`gap!(0N;0N;0n)];
        res:@[.commodity.calibrateCurve[;`schwartz2;rvCalCfgL];select tenor,price from sub;{[e] `rvError}];
        if[-11h=type res; :`longYM`shortYM`gap!(0N;0N;0n)];
        err:(res`perTenorError)`error;
        contracts:sub`contractYM;
        `longYM`shortYM`gap!(contracts err?max err;contracts err?min err;(max err)-min err)
        }[signalCurveHist;rvCalCfg;];
    rvList:rvIdent each dates;
    rvLongSeq:rvList[;`longYM]; rvShortSeq:rvList[;`shortYM];
    rvSignal:0f^rvList[;`gap];
    rvSpreadReturn:0n,{[priceMap;dates;rvLongSeq;rvShortSeq;t]
        ly:rvLongSeq t-1; sy:rvShortSeq t-1;
        lpNow:priceMap[dates t]ly; lpPrev:priceMap[dates t-1]ly;
        spNow:priceMap[dates t]sy; spPrev:priceMap[dates t-1]sy;
        $[(null ly)|(null sy)|(null lpNow)|(null lpPrev)|(null spNow)|null spPrev; 0n;
          ((lpNow%lpPrev)-1f)-((spNow%spPrev)-1f)]}[priceMap;dates;rvLongSeq;rvShortSeq;] each 1+til nDates-1;
    volTargetScaleRV:volScaleOf 0f^rvSpreadReturn;
    / --- assemble standard-schema path + signal columns ---
    pathTbl:flip `stepIndex`stepDate`spot`volatility`riskFreeRate`dividendYield`marketPrice`status`frontReturn`frontPrice`frontContractYM`daysToExpiry`curveSlopeCarry`momentum`chi`xi`chiZ`xiMomentum`convenienceYield`deferredReturn`nearFarSpreadReturn`rvSignal`rvSpreadReturn`isTrain`volTargetScale`volTargetScaleNearFar`volTargetScaleRV!(
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
        xiMomentum;
        convenienceYield;
        0f^deferredReturn;
        nearFarSpreadReturn;
        rvSignal;
        0f^rvSpreadReturn;
        isTrain;
        nDates#volTargetScale;
        nDates#volTargetScaleNearFar;
        nDates#volTargetScaleRV);
    `path`kalmanParams`trainEndDate`nTrain`nTest`volTargetScale`volTargetScaleNearFar`volTargetScaleRV`stationaryStd`heldContracts`productTag`deseasonalize`monthFactors!(
        pathTbl;kalParams;trainEndDate;sum isTrain;sum not isTrain;volTargetScale;volTargetScaleNearFar;volTargetScaleRV;stationaryStd;distinct heldSeq;cfg`productTag;deseasonalize;monthFactors)
 };

/ ============================================================================
/ Shared return-backtest core (vol-targeted position * return P&L, futures MTM)
/ ============================================================================
/ Strategies map a signal to a raw target in {-1,0,+1} (or continuous); the core
/ scales by volTargetScale*notional, lags the position one step, charges turnover
/ cost, and realises P&L to cash. PV=cash so deltaPV==stepPnl. Strategies differ
/ only in the raw-target rule; the accounting is shared and identical.

.strategy.commodityBT.__rowEmitCols:`stepIndex`stepDate`frontPrice`signal`rawTarget`position`frontReturn`positionPnl`turnoverCost`stepPnl`cumulativePnl`isTrain`status`message`targetPosition`proportionalCost`slippageCost`fixedCost`filledTurnover`targetTurnover;

/ The traded return / vol-scale columns default to the front series; spread
/ strategies (B2/B4) override them via returnColumn / volScaleColumn in stratCfg.
.strategy.commodityBT.__returnColumn:{[stratCfg] $[`returnColumn in key stratCfg; stratCfg`returnColumn; `frontReturn]};
.strategy.commodityBT.__volScaleColumn:{[stratCfg] $[`volScaleColumn in key stratCfg; stratCfg`volScaleColumn; `volTargetScale]};

/ Build the execution ctx for a step (8-param lambda cap -> bundle in a dict).
.strategy.commodityBT.__execCtx:{[stepRow;prevPos]
    `refPrice`barVolume`volatility`currentPos!(
        stepRow`frontPrice;
        $[`volume in key stepRow; stepRow`volume; 0n];
        stepRow`volatility;
        prevPos)
 };

.strategy.commodityBT.coreInit:{[trade;firstStep;stratCfg;rawTarget0;signalVal0]
    notional:trade`notional;
    scale:firstStep .strategy.commodityBT.__volScaleColumn stratCfg;
    position0:rawTarget0*scale*notional;
    execCfg:.exec.__resolve stratCfg;
    fillRes:.exec.fill[position0;.strategy.commodityBT.__execCtx[firstStep;0f];execCfg];
    filled0:fillRes`filledQty;
    entryCost:fillRes`totalCost;
    stepPnl0:neg entryCost;
    rowEmit:.strategy.commodityBT.__rowEmitCols!(
        firstStep`stepIndex;firstStep`stepDate;firstStep`frontPrice;signalVal0;rawTarget0;filled0;
        firstStep .strategy.commodityBT.__returnColumn stratCfg;0f;entryCost;stepPnl0;stepPnl0;firstStep`isTrain;`OK;"";
        position0;fillRes`proportionalCost;fillRes`slippageCost;fillRes`fixedCost;abs filled0;abs position0);
    `cash`prevPosition`prevRawTarget`cumulativePnl`notional`rowEmit!(
        neg entryCost;filled0;rawTarget0;stepPnl0;notional;rowEmit)
 };

.strategy.commodityBT.coreStep:{[state;marketStep;stratCfg;rawTarget;signalVal]
    notional:state`notional;
    scale:marketStep .strategy.commodityBT.__volScaleColumn stratCfg;
    prevPos:state`prevPosition;
    ret:marketStep .strategy.commodityBT.__returnColumn stratCfg;
    desired:rawTarget*scale*notional;
    positionPnl:prevPos*ret;
    order:desired-prevPos;
    execCfg:.exec.__resolve stratCfg;
    fillRes:.exec.fill[order;.strategy.commodityBT.__execCtx[marketStep;prevPos];execCfg];
    filledQty:fillRes`filledQty;
    newPos:prevPos+filledQty;
    totalCost:fillRes`totalCost;
    stepPnl:positionPnl-totalCost;
    newCash:(state`cash)+stepPnl;
    cumulativePnl:(state`cumulativePnl)+stepPnl;
    rowEmit:.strategy.commodityBT.__rowEmitCols!(
        marketStep`stepIndex;marketStep`stepDate;marketStep`frontPrice;signalVal;rawTarget;newPos;
        ret;positionPnl;totalCost;stepPnl;cumulativePnl;marketStep`isTrain;`OK;"";
        desired;fillRes`proportionalCost;fillRes`slippageCost;fillRes`fixedCost;abs filledQty;abs order);
    @[state;`cash`prevPosition`prevRawTarget`cumulativePnl`rowEmit;:;(newCash;newPos;rawTarget;cumulativePnl;rowEmit)]
 };

/ Time-series performance summary (split into in-sample / out-of-sample by isTrain).
.strategy.commodityBT.coreSummary:{[resultTable;stratCfg;strategyName]
    annDays:$[`annualizationDays in key stratCfg; stratCfg`annualizationDays; 252f];
    notional:$[`notional in key stratCfg; stratCfg`notional; 1f];
    base:`strategyName`steps`totalPnl`testPnl`testAnnualReturn`testAnnualVol`testSharpe`grossSharpe`testMaxDrawdown`testHitRate`trainSharpe`turnoverCostTotal`proportionalCostTotal`slippageCostTotal`fixedCostTotal`filledTurnoverTotal`targetTurnoverTotal`fillRatio`status`errorMessage!(
        strategyName;0;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;`ERROR;"empty");
    if[(0=count resultTable)|not 98h=type resultTable; :base];
    okRows:resultTable where (resultTable`status)=`OK;
    if[0=count okRows; :@[base;`steps;:;count resultTable]];
    perf:.strategy.commodityBT.__perf[;annDays;notional];
    testRows:okRows where not okRows`isTrain;
    trainRows:okRows where okRows`isTrain;
    testPerf:perf testRows`stepPnl;
    trainPerf:perf trainRows`stepPnl;
    / GROSS = before execution cost: stepPnl + the booked totalCost (turnoverCost).
    grossPerf:perf (testRows`stepPnl)+testRows`turnoverCost;
    filledTot:sum okRows`filledTurnover;
    targetTot:sum okRows`targetTurnover;
    `strategyName`steps`totalPnl`testPnl`testAnnualReturn`testAnnualVol`testSharpe`grossSharpe`testMaxDrawdown`testHitRate`trainSharpe`turnoverCostTotal`proportionalCostTotal`slippageCostTotal`fixedCostTotal`filledTurnoverTotal`targetTurnoverTotal`fillRatio`status`errorMessage!(
        strategyName;count okRows;sum okRows`stepPnl;testPerf`totalPnl;testPerf`annualReturn;testPerf`annualVol;
        testPerf`sharpe;grossPerf`sharpe;testPerf`maxDrawdown;testPerf`hitRate;trainPerf`sharpe;sum okRows`turnoverCost;
        sum okRows`proportionalCost;sum okRows`slippageCost;sum okRows`fixedCost;filledTot;targetTot;
        $[targetTot>0f;filledTot%targetTot;1f];`OK;"")
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
    .cfg.strategy.convenienceYieldCarry
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
    .cfg.strategy.chiReversion
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
/ B1. timeSeriesMomentum - the commodity trend premium (sign of trailing return)
/ ============================================================================
.strategy.timeSeriesMomentum.defaultConfig:{[]
    .cfg.strategy.timeSeriesMomentum
 };
.strategy.timeSeriesMomentum.__rawTarget:{[mom;stratCfg]
    m:stratCfg`momentumMargin; $[mom>m;1f;mom<neg m;-1f;0f]
 };
.strategy.timeSeriesMomentum.init:{[trade;firstStep;model;fdmConfig;stratCfg]
    mom:firstStep`momentum;
    .strategy.commodityBT.coreInit[trade;firstStep;stratCfg;.strategy.timeSeriesMomentum.__rawTarget[mom;stratCfg];mom]
 };
.strategy.timeSeriesMomentum.step:{[state;marketStep;trade;model;fdmConfig;stratCfg]
    mom:marketStep`momentum;
    .strategy.commodityBT.coreStep[state;marketStep;stratCfg;.strategy.timeSeriesMomentum.__rawTarget[mom;stratCfg];mom]
 };
.strategy.timeSeriesMomentum.summary:{[resultTable;stratCfg]
    .strategy.commodityBT.coreSummary[resultTable;stratCfg;`timeSeriesMomentum]
 };
.strategy.register[`timeSeriesMomentum;.strategy.timeSeriesMomentum.init;.strategy.timeSeriesMomentum.step;.strategy.timeSeriesMomentum.summary;.strategy.timeSeriesMomentum.defaultConfig];

/ ============================================================================
/ B3. twoTimescale - trend-follow the slow factor xi + mean-revert the fast chi
/ ============================================================================
/ Combines a thresholded chi-reversion signal (fast timescale) with xi-momentum
/ trend-following (slow timescale), config-weighted. (Uses a thresholded rather
/ than hysteretic chi signal so the combined target is stateless.)
.strategy.twoTimescale.defaultConfig:{[]
    .cfg.strategy.twoTimescale
 };
.strategy.twoTimescale.__rawTarget:{[chiZ;xiMom;stratCfg]
    entryZ:stratCfg`entryZ;
    revertRaw:$[chiZ<neg entryZ;1f;chiZ>entryZ;-1f;0f];
    trendRaw:signum xiMom;
    ((stratCfg`revertWeight)*revertRaw)+(stratCfg`trendWeight)*trendRaw
 };
.strategy.twoTimescale.init:{[trade;firstStep;model;fdmConfig;stratCfg]
    raw:.strategy.twoTimescale.__rawTarget[firstStep`chiZ;firstStep`xiMomentum;stratCfg];
    .strategy.commodityBT.coreInit[trade;firstStep;stratCfg;raw;firstStep`chiZ]
 };
.strategy.twoTimescale.step:{[state;marketStep;trade;model;fdmConfig;stratCfg]
    raw:.strategy.twoTimescale.__rawTarget[marketStep`chiZ;marketStep`xiMomentum;stratCfg];
    .strategy.commodityBT.coreStep[state;marketStep;stratCfg;raw;marketStep`chiZ]
 };
.strategy.twoTimescale.summary:{[resultTable;stratCfg]
    .strategy.commodityBT.coreSummary[resultTable;stratCfg;`twoTimescale]
 };
.strategy.register[`twoTimescale;.strategy.twoTimescale.init;.strategy.twoTimescale.step;.strategy.twoTimescale.summary;.strategy.twoTimescale.defaultConfig];

/ ============================================================================
/ B4. storageCashCarry - harvest contango above storage cost via a calendar spread
/ ============================================================================
/ When contango (negative curve carry) exceeds the storage cost, put on a long-
/ near / short-far calendar spread (trade the precomputed nearFarSpreadReturn);
/ flat otherwise. The spread is spot-neutral, so it isolates the carry.
.strategy.storageCashCarry.defaultConfig:{[]
    .cfg.strategy.storageCashCarry
 };
.strategy.storageCashCarry.__rawTarget:{[curveCarry;stratCfg]
    contango:neg curveCarry;
    $[contango>(stratCfg`storageCostRate)+stratCfg`carryMargin;1f;0f]
 };
.strategy.storageCashCarry.init:{[trade;firstStep;model;fdmConfig;stratCfg]
    cc:firstStep`curveSlopeCarry;
    .strategy.commodityBT.coreInit[trade;firstStep;stratCfg;.strategy.storageCashCarry.__rawTarget[cc;stratCfg];cc]
 };
.strategy.storageCashCarry.step:{[state;marketStep;trade;model;fdmConfig;stratCfg]
    cc:marketStep`curveSlopeCarry;
    .strategy.commodityBT.coreStep[state;marketStep;stratCfg;.strategy.storageCashCarry.__rawTarget[cc;stratCfg];cc]
 };
.strategy.storageCashCarry.summary:{[resultTable;stratCfg]
    .strategy.commodityBT.coreSummary[resultTable;stratCfg;`storageCashCarry]
 };
.strategy.register[`storageCashCarry;.strategy.storageCashCarry.init;.strategy.storageCashCarry.step;.strategy.storageCashCarry.summary;.strategy.storageCashCarry.defaultConfig];

/ ============================================================================
/ B5. carryMomentumCombo - the documented backbone (carry + momentum, weighted)
/ ============================================================================
.strategy.carryMomentumCombo.defaultConfig:{[]
    .cfg.strategy.carryMomentumCombo
 };
.strategy.carryMomentumCombo.__rawTarget:{[cy;mom;stratCfg]
    rate:stratCfg`riskFreeRate; cm:stratCfg`carryMargin; mm:stratCfg`momentumMargin;
    carryRaw:$[cy>rate+cm;1f;cy<rate-cm;-1f;0f];
    momRaw:$[mom>mm;1f;mom<neg mm;-1f;0f];
    ((stratCfg`carryWeight)*carryRaw)+(stratCfg`momentumWeight)*momRaw
 };
.strategy.carryMomentumCombo.init:{[trade;firstStep;model;fdmConfig;stratCfg]
    raw:.strategy.carryMomentumCombo.__rawTarget[firstStep`convenienceYield;firstStep`momentum;stratCfg];
    .strategy.commodityBT.coreInit[trade;firstStep;stratCfg;raw;firstStep`convenienceYield]
 };
.strategy.carryMomentumCombo.step:{[state;marketStep;trade;model;fdmConfig;stratCfg]
    raw:.strategy.carryMomentumCombo.__rawTarget[marketStep`convenienceYield;marketStep`momentum;stratCfg];
    .strategy.commodityBT.coreStep[state;marketStep;stratCfg;raw;marketStep`convenienceYield]
 };
.strategy.carryMomentumCombo.summary:{[resultTable;stratCfg]
    .strategy.commodityBT.coreSummary[resultTable;stratCfg;`carryMomentumCombo]
 };
.strategy.register[`carryMomentumCombo;.strategy.carryMomentumCombo.init;.strategy.carryMomentumCombo.step;.strategy.carryMomentumCombo.summary;.strategy.carryMomentumCombo.defaultConfig];

/ ============================================================================
/ B2. curveRelativeValue - long the cheapest / short the richest tenor (by the
/ schwartz2 calibration residual), expecting convergence. Trades the precomputed
/ rvSpreadReturn (cheap-leg return minus rich-leg return), level-neutral.
/ ============================================================================
.strategy.curveRelativeValue.defaultConfig:{[]
    .cfg.strategy.curveRelativeValue
 };
.strategy.curveRelativeValue.__rawTarget:{[gap;stratCfg] $[gap>stratCfg`minGap;1f;0f]};
.strategy.curveRelativeValue.init:{[trade;firstStep;model;fdmConfig;stratCfg]
    g:firstStep`rvSignal;
    .strategy.commodityBT.coreInit[trade;firstStep;stratCfg;.strategy.curveRelativeValue.__rawTarget[g;stratCfg];g]
 };
.strategy.curveRelativeValue.step:{[state;marketStep;trade;model;fdmConfig;stratCfg]
    g:marketStep`rvSignal;
    .strategy.commodityBT.coreStep[state;marketStep;stratCfg;.strategy.curveRelativeValue.__rawTarget[g;stratCfg];g]
 };
.strategy.curveRelativeValue.summary:{[resultTable;stratCfg]
    .strategy.commodityBT.coreSummary[resultTable;stratCfg;`curveRelativeValue]
 };
.strategy.register[`curveRelativeValue;.strategy.curveRelativeValue.init;.strategy.curveRelativeValue.step;.strategy.curveRelativeValue.summary;.strategy.curveRelativeValue.defaultConfig];

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

/ ============================================================================
/ Part A (v0.54). Walk-forward robustness: the OOS verdict as a DISTRIBUTION.
/ ----------------------------------------------------------------------------
/ A single year's Sharpe has a standard error on the order of +/-1, so a point
/ estimate is not a verdict. walkForward runs the existing per-run backtest over
/ multiple sequential train/test splits (rolling or expanding), each estimating
/ the Kalman params on THAT split's train only (causal, no cross-split leak),
/ and aggregates each strategy's OOS metrics into mean +/- dispersion + how many
/ splits were positive. Splits are INDEX-BASED over the distinct curve dates, so
/ appending future data adds new splits without changing existing ones (the
/ cross-split causality guarantee). Spans are in number-of-dates.
/ ============================================================================

.strategy.commodityBT.defaultSplitCfg:{[]
    `scheme`trainSpan`testSpan`maxSplits!(`expanding;120;60;8)
 };

/ Index-based split specs over nDates dated points. Returns a list of dicts
/ (trainStartIdx, trainEndIdx, testEndIdx); stops when a test window runs past
/ the data. Expanding: train grows from 0; rolling: fixed-width train slides.
.strategy.commodityBT.__splits:{[scheme;nDates;trainSpan;testSpan;maxSplits]
    splits:();
    i:0;
    while[i<maxSplits;
        trainStartIdx:$[scheme=`rolling; i*testSpan; 0];
        trainEndIdx:$[scheme=`rolling; trainStartIdx+trainSpan-1; trainSpan+(i*testSpan)-1];
        testEndIdx:trainEndIdx+testSpan;
        if[testEndIdx>=nDates; :splits];
        splits,:enlist `trainStartIdx`trainEndIdx`testEndIdx!(trainStartIdx;trainEndIdx;testEndIdx);
        i+:1];
    splits
 };

/ Aggregate the per-split detail (one row per strategy per split) into a per-
/ strategy distribution. avg / dev skip nulls (strategies that never traded).
.strategy.commodityBT.__aggregateSplits:{[detail]
    if[(0=count detail)|not 98h=type detail;
        :([] strategyName:0#`; nSplits:0#0; nTraded:0#0; meanSharpe:0#0Nf; sharpeStd:0#0Nf; splitsPositive:0#0; meanAnnReturn:0#0Nf; meanMaxDrawdown:0#0Nf)];
    0!select
        nSplits:count i,
        nTraded:sum not null testSharpe,
        meanSharpe:avg testSharpe,
        sharpeStd:dev testSharpe,
        splitsPositive:sum testSharpe>0f,
        meanAnnReturn:avg testAnnualReturn,
        meanMaxDrawdown:avg testMaxDrawdown
        by strategyName from detail
 };

.strategy.commodityBT.walkForward:{[strategies;trade;curveHistory;sigCfg;splitCfg]
    if[not 98h=type curveHistory; '"walkForward: curveHistory must be a table"];
    cfg:.strategy.commodityBT.defaultSplitCfg[];
    if[count splitCfg; cfg:cfg,splitCfg];
    dates:asc distinct curveHistory`asofDate;
    nDates:count dates;
    splits:.strategy.commodityBT.__splits[cfg`scheme;nDates;cfg`trainSpan;cfg`testSpan;cfg`maxSplits];
    if[0=count splits; '"walkForward: no valid splits for the given spans and history window"];
    / Bundle fixed args into a ctx dict (8-param lambda cap).
    ctx:`strategies`trade`curveHistory`sigCfg`dates`model`fdmCfg!(
        strategies;trade;curveHistory;sigCfg;dates;.model.createBlackScholesModel[];()!());
    runSplit:{[spec;ctx;splitId]
        dts:ctx`dates;
        ts:dts spec`trainStartIdx; te:dts spec`trainEndIdx; xe:dts spec`testEndIdx;
        subHist:select from (ctx`curveHistory) where asofDate within (ts;xe);
        sig:.strategy.path.commoditySignals[subHist;@[ctx`sigCfg;`trainEndDate;:;te]];
        suite:.strategy.commodityBT.runSuite[ctx`strategies;ctx`trade;sig`path;ctx`model;ctx`fdmCfg;()!()];
        perf:select strategyName,testSharpe,testAnnualReturn,testMaxDrawdown,testHitRate from suite`ranked;
        update splitId:splitId, trainEnd:te, testEnd:xe from perf};
    detail:raze runSplit[;ctx;]'[splits;til count splits];
    splitDates:([] splitId:til count splits;
        trainStart:dates splits[;`trainStartIdx]; trainEnd:dates splits[;`trainEndIdx]; testEnd:dates splits[;`testEndIdx]);
    `aggregate`detail`splits!(.strategy.commodityBT.__aggregateSplits detail;detail;splitDates)
 };

/ ============================================================================
/ Part C (v0.55). Cross-commodity robustness: aggregate walk-forward OOS across
/ BOTH commodities and all splits into one robustness score per strategy.
/ ----------------------------------------------------------------------------
/ The verdict is which strategies are robustly positive across commodities AND
/ windows (a real claim to edge) vs commodity/period-specific (noise). Ranked by
/ mean OOS Sharpe across all (commodity x split) cells; fractionPositive is the
/ share of cells with positive Sharpe.
/ ============================================================================
.strategy.commodityBT.crossCommodity:{[detailByCommodity]
    if[0=count detailByCommodity; '"crossCommodity: empty detailByCommodity"];
    commodities:key detailByCommodity;
    tagged:raze {[d;c] update commodity:c from d}'[value detailByCommodity;commodities];
    agg:0!select
        nCells:count i,
        nTraded:sum not null testSharpe,
        meanSharpe:avg testSharpe,
        sharpeStd:dev testSharpe,
        fractionPositive:(sum testSharpe>0f)%count i,
        meanAnnReturn:avg testAnnualReturn
        by strategyName from tagged;
    `meanSharpe xdesc agg
 };
