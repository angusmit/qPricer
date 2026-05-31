/ backtest/replay.q - the event-driven replay engine (.backtest.replay.*) (v0.77, Research OS R12)
/ ----------------------------------------------------------------------------
/ ARCHITECTURE.md Part II 11.8(d) / 13-R12. THE KEYSTONE of the evidence layer. R9 gave the single
/ as-of door, R10 the curve, R11 roll discipline - all as-of-clean inputs. R12 folds them into a
/ DETERMINISTIC date-by-date loop that produces REALISTIC, cost-and-roll-aware PnL on the ACTUAL
/ contracts, plus the auditable run record R13's evidence audit will bite on. After R12, replay mode
/ is the promotion gate: every verdict finally corresponds to something a desk could have traded.
/ -
/ THE DESIGN (bounded + faithful). The replay REUSES the causal signal->target (the unchanged
/ .strategy.path.commoditySignals + the strategy's research-mode run produce the per-day target
/ position causally - the R6 faithful-delegation precedent) and INDEPENDENTLY RE-DERIVES the as-of
/ execution: the active contract (R11's days_before_expiry engine, parameterised to the strategy's
/ own roll convention with W=1 so it tracks the research heldSeq), the as-of return (R9's door - the
/ active contract's price at d vs the prior date, both date<=d), the fills (the UNCHANGED .exec.fill
/ + opt-in realism), the position book, the PnL (price + cost + financing on the actual contract),
/ and the run record. So it detects EXECUTION / PRICE / ROLL look-ahead (the R12 concern) and
/ reproduces the research-mode PnL to a tight tolerance (the look-ahead detector, Part 4).
/ -
/ BYTE-IDENTITY: this is a NEW path. It does NOT modify .exec.fill or the research-mode backtester;
/ the realism levers (.cfg.replay) are all OFF by default, so a frictionless replay calls the
/ UNCHANGED .exec.fill and the canonicals stay byte-identical. ENGINE INFRASTRUCTURE (like the
/ research-mode backtester) - NOT registered as an R2 capability and NOT carded.
/ -
/ LAYER: HIGH within backtest/ - composes state/ (R9), curve/ (R10), roll/ (R11), execution/, and
/ the signal/strategy path, ALL below or beside it (downward only). Loads after commodityStrategies.q;
/ opens nothing at import.
/ -
/ RESERVED-NAME NOTE: parameter `asOf` NEVER `asof`; `comm` NOT `commodity`. Builtins not shadowed
/ (sum/sums/deltas/next/prev/first/last/over/scan/get/key/value/cols).
/ ----------------------------------------------------------------------------

.backtest.replay.defaultConfig:{[] .cfg.replay};

/ ── inputs: the causal path + the research-mode run (reused faithfully) ───────

/ Build the shared inputs once: the causal signal path, the strategy's research-mode result (whose
/ per-day targetPosition is the causal signal->target the replay re-executes), the trade + configs.
/ Both .backtest.replay.run and .backtest.replay.faithfulness use this, so they see identical inputs.
.backtest.replay.__buildInputs:{[strategyName;comm;userCfg]
    curveHist:.data.hdb.curveHistory[comm; .data.hdb.dates comm];
    sigCfg:`rollDaysBeforeExpiry`momentumLookback`txnCostRate`targetVol`riskFreeRate`kalmanEstCfg`annualizationDays!(
        5;20;0.0005;0.15;0.02;(`gridSteps`refineRounds`nSweeps!(7;3;3));252f);
    if[`sigCfg in key userCfg; sigCfg:sigCfg,userCfg`sigCfg];
    path:(.strategy.path.commoditySignals[curveHist;sigCfg])`path;
    trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
        `RPLAY;`WTI;`equityOption;`european;`call;60f;0.25;1f);
    stratCfg:.strategy.defaultConfig strategyName;
    if[`stratCfg in key userCfg; stratCfg:stratCfg,userCfg`stratCfg];
    res:.strategy.runAndSummarize[strategyName;trade;path;(.model.createBlackScholesModel[]);(()!());stratCfg];
    `path`res`trade`stratCfg`sigCfg!(path;res;trade;stratCfg;sigCfg)
 };

/ ── the deterministic fold: a pure step folded over the trading dates ─────────

/ One step: given the running state and one date-row, build the as-of slice (R9), the active
/ contract (R11), re-derive the as-of return, generate the order (target - current + carry), fill it
/ (the unchanged .exec.fill + opt-in realism), book PnL (price - cost - financing), append the record.
/ Pure given (state; row) + the HDB - reproducible + fully auditable (the whole trail is the record).
.backtest.replay.__step:{[state;row]
    comm:state`comm; replayCfg:state`replayCfg; rollCfg:state`rollCfg;
    d:row`stepDate; prevDate:row`prevDate;
    pos:state`pos; cash:state`cash; carry:state`carry; cumPnl:state`cumPnl; activePrev:state`activePrev;
    / R9 door: ONLY data knowable as of d, plus the provenance stamp (the audit cross-reference).
    a:.state.asof[d;comm]; slice:a`data; prov:a`provenance;
    / The active contract is decided as of the PRIOR trading date (refDate), then traded on d - the
    / causally-honest convention (you know which contract you hold before the day opens), and the one
    / the research path uses: its heldSeq is a one-day-lagged roll scan, so deciding the roll on refDate
    / makes the replay's active track heldSeq exactly (faithfulness). refDate<=d, so it stays as-of-safe.
    / The live universe is the SAME computation as .roll.__liveAsOf (reusing the slice we already pulled
    / through the door), with a LARGE volume window so the days_before_expiry rule sees EVERY live
    / contract by expiry; R11's rule engine then picks the active contract (blend W=1).
    refDate:$[null prevDate; d; prevDate];
    effDate:max slice`date;
    live:`expiry xasc 0!select expiry:first expiry, vol:sum volume, effDate:max date
        by contractYM from select from slice where expiry>refDate, date within (effDate-state`volWin;effDate);
    c:$[0=count live; activePrev; (.roll.__daysBeforeExpiry[refDate;live;rollCfg])`active];
    / as-of return on the ACTUAL active contract (price at d vs the prior date - both date<=d).
    priceNow:first exec settle from slice where contractYM=c, date=d;
    pricePrev:$[null prevDate; 0Nf; first exec settle from slice where contractYM=c, date=prevDate];
    retAsOf:$[(null prevDate) or (null priceNow) or null pricePrev; 0f; (priceNow%pricePrev)-1f];
    barVol:first exec volume from slice where contractYM=c, date=d;
    isRoll:(not null activePrev) and not c=activePrev;
    / order = causal target - current position + any carried (unfilled) remainder.
    desired:row`targetPosition;
    order:(desired-pos)+carry;
    / opt-in participation cap (default off -> filled=order, remainder=0).
    capRes:.exec.participationCap[order;barVol;replayCfg`participationRate];
    filledOrder:capRes`filled; rem:capRes`remainder;
    ctx:`refPrice`barVolume`volatility`currentPos!(priceNow;barVol;0Nf;pos);
    fillRes:.exec.fill[filledOrder;ctx;state`execCfgBase];   / the UNCHANGED chokepoint
    filledQty:fillRes`filledQty;
    newPos:pos+filledQty;
    positionPnl:pos*retAsOf;
    / opt-in roll cost (the held exposure rolled on a roll date) + financing (default both 0).
    rollCost:.exec.rollPenalty[pos*$[null priceNow;0f;priceNow]; isRoll; replayCfg`rollPenaltyBps];
    financingCost:((replayCfg`financingRate)%replayCfg`annualizationDays)*(abs pos)*$[null priceNow;0f;priceNow];
    totalCost:(fillRes`totalCost)+rollCost;
    stepPnl:positionPnl-totalCost-financingCost;
    newCash:cash+stepPnl; newCum:cumPnl+stepPnl;
    newCarry:$[(replayCfg`remainderPolicy)=`carry; rem; 0f];
    rec:`stepIndex`stepDate`asOf`activeContract`isRoll`rollFrom`rawTarget`targetPosition`order`filledQty`fillPrice`refPrice`frontReturn`barVolume`proportionalCost`slippageCost`fixedCost`rollCost`financingCost`totalCost`positionPnl`stepPnl`cumPnl`position`cash`provNRows`provDateTo!(
        row`stepIndex; d; d; c; isRoll; $[isRoll;activePrev;0Nj]; row`rawTarget; desired; order; filledQty;
        fillRes`fillPrice; priceNow; retAsOf; barVol; fillRes`proportionalCost; fillRes`slippageCost;
        fillRes`fixedCost; rollCost; financingCost; totalCost; positionPnl; stepPnl; newCum; newPos; newCash;
        prov`nRows; prov`dateTo);
    @[state; `cash`pos`carry`cumPnl`activePrev`rowEmit; :; (newCash;newPos;newCarry;newCum;c;rec)]
 };

/ The fold core: fold .backtest.replay.__step over the sorted target rows and assemble the auditable
/ run record (`meta`steps`rollEvents`provenance). rowsIn carries per-date stepIndex/stepDate/rawTarget/
/ targetPosition (the causal signal->target). foldCfg bundles replayCfg/rollDays/execCfgBase/notional/
/ fromDate/toDate/label/commodity (bundled to stay within the 8-param limit). The active contract +
/ as-of return + fills + PnL are RE-DERIVED inside __step from the as-of data. fromDate/toDate default
/ (0Nd) to the full range. Deterministic (no RNG); two runs are identical. (The test drives this with
/ synthetic rows + a synthetic `futures` global, bypassing the heavy signal builder.)
.backtest.replay.__fold:{[rowsIn;foldCfg]
    comm:foldCfg`commodity; replayCfg:foldCfg`replayCfg;
    rollCfg:`type`rollDays`W`priceSource`method!(`days_before_expiry; foldCfg`rollDays; replayCfg`rollW; `settle; `difference);
    rows:update prevDate:prev stepDate from `stepDate xasc rowsIn;
    fromD:$[null foldCfg`fromDate; min rows`stepDate; foldCfg`fromDate];
    toD:$[null foldCfg`toDate; max rows`stepDate; foldCfg`toDate];
    rows:select from rows where stepDate within (fromD;toD);
    if[0=count rows; '"replay.__fold: no dates in the requested range"];
    seed:`cash`pos`carry`cumPnl`activePrev`comm`replayCfg`execCfgBase`rollCfg`volWin`notional`rowEmit!(
        0f; 0f; 0f; 0f; 0Nj; comm; replayCfg; foldCfg`execCfgBase; rollCfg; 1000000; foldCfg`notional; ()!());
    states:.backtest.replay.__step\[seed; 0!rows];
    steps:.strategy.__rowDictsToTable states[;`rowEmit];
    rollEvents:select date:stepDate, commodity:comm, fromContract:rollFrom, toContract:activeContract, reason:`replayRoll from steps where isRoll;
    provenance:select stepDate, asOf, provNRows, provDateTo from steps;
    metaDict:`strategy`commodity`fromDate`toDate`nSteps`totalPnl`grossPnl`totalCost`endPosition`endCash`nRolls!(
        foldCfg`label; comm; fromD; toD; count steps; sum steps`stepPnl; (sum steps`stepPnl)+sum steps`totalCost;
        sum steps`totalCost; last steps`position; last steps`cash; count rollEvents);
    `meta`steps`rollEvents`provenance!(metaDict;steps;rollEvents;provenance)
 };

/ Fold from already-built inputs (so faithfulness can build the inputs ONCE and reuse them).
.backtest.replay.__runFromInputs:{[inp;comm;replayCfg;fromDate;toDate;strategyName]
    okRes:select from (inp`res)`result where status=`OK;
    if[0=count okRes; '"replay.run: no OK research rows for ",string strategyName];
    rows:select stepIndex, stepDate, rawTarget, targetPosition from okRes;
    foldCfg:`commodity`replayCfg`rollDays`execCfgBase`notional`fromDate`toDate`label!(
        comm; replayCfg; (inp`sigCfg)`rollDaysBeforeExpiry; .exec.__resolve inp`stratCfg; (inp`trade)`notional; fromDate; toDate; strategyName);
    .backtest.replay.__fold[rows;foldCfg]
 };

/ The replay run: reuse the causal signal->target (the unchanged strategy research run), then fold the
/ event-driven, as-of-verified engine over the dates. fromDate/toDate default (0Nd) to the full range.
.backtest.replay.run:{[strategyName;comm;fromDate;toDate;userCfg]
    inp:.backtest.replay.__buildInputs[strategyName;comm;userCfg];
    replayCfg:.cfg.replay; if[`replay in key userCfg; replayCfg:replayCfg,userCfg`replay];
    .backtest.replay.__runFromInputs[inp;comm;replayCfg;fromDate;toDate;strategyName]
 };

/ ── faithfulness (the look-ahead detector) ───────────────────────────────────

/ Run the research mode and a FRICTIONLESS replay of the same strategy and compare their per-day PnL.
/ The comparison is over the dates where the replay's active contract AGREES with the research path's
/ held contract - i.e. where the roll map names the SAME contract. On those dates the frictionless
/ replay must reproduce the research PnL to a tight RELATIVE tolerance; a divergence there is look-ahead
/ or a roll/fill bug (the detector). It is a tolerance, NOT an equality, because the replay re-derives
/ the per-day return from the as-of slice (R9) INDEPENDENTLY of the path builder's vectorised
/ construction, so float rounding order can differ in the last bits. The few DATA-EDGE dates where the
/ contracts differ (the research path holds a contract that has stopped trading at the tail of the
/ data, while the replay rolls to the only listed contract) are reported SEPARATELY as `nEdge`/
/ `edgeDates - not hidden, not a look-ahead. The research path is UNTOUCHED (canonical totalPnl is
/ byte-identical). tol is the tight relative tolerance for the agree-dates total.
.backtest.replay.faithfulness:{[strategyName;comm;tol]
    inp:.backtest.replay.__buildInputs[strategyName;comm;()!()];
    rr:.backtest.replay.__runFromInputs[inp;comm;.cfg.replay;0Nd;0Nd;strategyName];
    res:select stepDate, rPnl:stepPnl from (inp`res)`result where status=`OK;
    held:select stepDate, heldYM:frontContractYM from inp`path;
    qy:select stepDate, qPnl:stepPnl, qActive:activeContract from rr`steps;
    cmp:0!((`stepDate xkey res) lj (`stepDate xkey qy) lj `stepDate xkey held);
    clean:select from cmp where qActive=heldYM;
    edge:select from cmp where qActive<>heldYM;
    researchTotal:sum clean`rPnl; replayTotal:sum clean`qPnl;
    relDiff:(abs replayTotal-researchTotal)%0f|1e-12|abs researchTotal;
    `researchTotal`replayTotal`relDiff`tol`pass`nEdge`edgeDates!(
        researchTotal; replayTotal; relDiff; tol; relDiff<tol; count edge; edge`stepDate)
 };

/ ── thin opt-in persistence (splayed, UNPARTITIONED replayRuns; §3 storage) ──

/ Persist a run's `steps to the splayed, UNPARTITIONED replayRuns (p# on commodity, sym columns via
/ .Q.en, written with set - the SAME storage decision as the HDB; NOT date-partitioned). A runId tags
/ the rows. gitignored; tests use the in-memory record. Returns the row count written.
.backtest.replay.persist:{[hdbPath;runId;runRecord]
    steps:update runId:runId, commodity:(runRecord[`meta]`commodity) from runRecord`steps;
    t:.Q.en[hsym `$hdbPath; update `p#commodity from `commodity`stepDate xasc steps];
    (hsym `$hdbPath,"/replayRuns/") set t;
    `rows!enlist count t
 };
.backtest.replay.open:{[hdbPath]
    sp:hsym `$hdbPath,"/replayRuns/";
    if[0=count key sp; '"replay.open: no replayRuns at ",hdbPath];
    get sp
 };

-1 "replay.q loaded - .backtest.replay.* event-driven replay engine ready (frictionless default; faithful to research mode)";
