/ templates/seasonal_calendar_spread.q - the seasonalCalendarSpread template (v0.81, Research OS R16)
/ ----------------------------------------------------------------------------
/ THE CAPSTONE. The first real research output: fade the crude calendar spread against its
/ SEASONALLY-ADJUSTED same-calendar-month z-score (R14) - SHORT the spread when it is unusually wide
/ for the season (z>+entry), LONG when unusually narrow (z<-entry), EXIT on a z-cross of zero or maxHold.
/ Trades the ACTUAL contracts (the front + deferred legs the curve names, R11/R10), books NET PnL via
/ the fill model (.exec.fill), and - unlike R6's relativeValue - emits an R12-shape RUN RECORD so the
/ evidence audit (R13) and the PnL attribution (R15) can consume it end-to-end. Edge: STRUCTURAL (the
/ seasonally-adjusted spread reverts to its seasonal norm). Composes R9/R10/R11/R14 downward.
/ -
/ THE POINT (read first): the foundation was built for an HONEST answer, not a passing one. The
/ parameters are PRE-REGISTERED in .cfg.strategy.seasonalCalSpread (FIXED, NOT swept) - tuning the legs
/ / z-window / thresholds is exactly where a calendar spread overfits, and the deflation + walk-forward
/ + sealed holdout exist to punish it. A "research, not tradeable" verdict is the system WORKING.
/ -
/ TEMPLATE-SPECIFIC GATE that BITES: the spread must MEAN-REVERT (reuse R6's .template.rv.stationarity,
/ AR(1)/OU). A random-walk spread -> fading the z is fading noise -> REJECT before the universal gates.
/ -
/ RESERVED-NAME NOTE: `asOf` NEVER `asof`; `comm` NOT `commodity`; calendar month via (`month$date) mod
/ 12 (.mm is atom-only); builtins not shadowed (dev/avg/var/first/last/signum/sums/deltas/prev).
/ ----------------------------------------------------------------------------

/ The calendar-spread series WITH the front/deferred leg contracts, from a curve history (asofDate/
/ tenor/price/contractYM): per date, front = nearest tenor, deferred = the deferredIdx-th (clamped).
.template.scs.__spreadSeries:{[ch;deferredIdx;rollBuffer]
    d2e:`int$(ch`expiry)-ch`asofDate;                            / days to expiry (mask, NOT qSQL on a local)
    ch:ch where d2e>rollBuffer;                                  / roll discipline: skip the expiring contract
    s:`asofDate`tenor xasc ch;
    g:0!select front:first price, deferred:price@deferredIdx&-1+count price,
        frontYM:first contractYM, deferredYM:contractYM@deferredIdx&-1+count contractYM by asofDate from s;
    `date`spread`front`frontYM`deferredYM!(g`asofDate; (g`front)-g`deferred; g`front; g`frontYM; g`deferredYM)
 };

/ The CAUSAL seasonal same-calendar-month z of the spread: z[t] vs the SAME-calendar-month spreads
/ STRICTLY BEFORE t (never the full sample). null below minObs same-month observations. PURE.
.template.scs.__seasonalZ:{[spread;dates;minObs]
    cm:(`month$dates) mod 12;
    {[spread;cm;minObs;i]
        hist:i#spread; histCm:i#cm;                          / strictly before t (causal)
        sm:hist where histCm=cm i;
        sd:$[1<count sm; dev sm; 0f];
        $[((count sm)>=minObs) and sd>0f; (spread[i]-avg sm)%sd; 0n]}[spread;cm;minObs;] each til count spread
 };

/ The position state machine over the z series: flat -> SHORT (-1) when z>+entry, LONG (+1) when
/ z<-entry; while in a trade, EXIT (->0) when z crosses zero (reverted) or maxHold days elapsed. A
/ scan threading (position; daysHeld). Null z (below min-N) is treated as no signal (stay/flat).
.template.scs.__positions:{[z;entry;maxHold]
    step:{[entry;maxHold;st;zi]
        pos:st 0; held:st 1;
        $[pos=0;
            $[null zi; (0;0);
              zi>entry; (-1;0);                              / spread too wide -> short it
              zi<neg entry; (1;0);                           / spread too narrow -> long it
              (0;0)];
            / in a trade: exit on z-cross-of-zero (reverted) or maxHold
            [reverted:(not null zi) and $[pos<0; zi<=0f; zi>=0f];
             $[reverted or held>=maxHold-1; (0;0); (pos;held+1)]]]}[entry;maxHold];
    (step\[(0;0); z])[;0]                                    / one state per z element (the seed is NOT in the scan output); take the position column
 };

/ The end-to-end as-of REPLAY of the spread: fold over the dates trading the legs, book NET PnL via
/ .exec.fill, emit an R12-shape run record (front leg as the nominal activeContract; spread DOLLAR PnL).
/ Returns `meta`steps`rollEvents`provenance (the shape R13's audit + R15's attribution consume).
.template.scs.__replay:{[comm;fromDate;toDate;cfg]
    ch0:.data.hdb.curveHistory[comm; .data.hdb.dates comm];
    ser:.template.scs.__spreadSeries[ch0; cfg`deferredIdx; cfg`rollBuffer];
    allDates:ser`date;
    fromD:$[null fromDate; min allDates; fromDate];          / null bounds -> the full range
    toD:$[null toDate; max allDates; toDate];
    dates:allDates;
    keep:dates within (fromD;toD);
    dates:dates where keep; spread:(ser`spread) where keep; front:(ser`front) where keep;
    frontYM:(ser`frontYM) where keep; deferredYM:(ser`deferredYM) where keep;
    z:.template.scs.__seasonalZ[spread;dates;cfg`minObs];
    pos:cfg[`notional]*.template.scs.__positions[z;cfg`entryThreshold;cfg`maxHold];
    n:count dates;
    posLag:0f,-1_pos;                                        / position held INTO each day's spread move
    dSpread:0f,1_deltas spread;                              / spread change (0 on day 0)
    positionPnl:posLag*dSpread;                              / spread dollar PnL
    order:deltas pos; filledQty:order;                       / spread position changes (full fill)
    execCfg:.exec.__resolve cfg;
    cost:{[execCfg;ord] (.exec.fill[ord; `refPrice`currentPos!(1f;0f); execCfg])`totalCost}[execCfg] each order;
    stepPnl:positionPnl-cost;
    isRoll:differ frontYM; isRoll[0]:0b;
    steps:flip `stepIndex`stepDate`asOf`activeContract`isRoll`rollFrom`position`frontReturn`positionPnl`proportionalCost`slippageCost`fixedCost`rollCost`financingCost`totalCost`stepPnl`cumPnl`filledQty`refPrice`provDateTo!(
        til n; dates; dates; frontYM; isRoll; ?[isRoll; prev frontYM; 0Nj]; pos; dSpread; positionPnl;
        cost; n#0f; n#0f; n#0f; n#0f; cost; stepPnl; sums stepPnl; filledQty; front; dates);
    rollEvents:select date:stepDate, commodity:comm, fromContract:rollFrom, toContract:activeContract, reason:`legRoll from steps where isRoll;
    `meta`steps`rollEvents`provenance!(
        `strategy`commodity`fromDate`toDate`nSteps`totalPnl`spreadAr1!(
            `seasonalCalSpread; comm; min dates; max dates; n; sum stepPnl; (.template.rv.stationarity spread)`ar1);
        steps; rollEvents; select stepDate, asOf, provDateTo from steps)
 };

/ The runner seam for the gov gates: run the template over [from;to] -> (date;pnl). Reuses __replay.
.template.scs.runner:{[comm;cfg] {[comm;cfg;from;to] r:.template.scs.__replay[comm;from;to;cfg]; select date:stepDate, pnl:stepPnl from r`steps}[comm;cfg]};

/ The template run (the `template` contract): inputs carry commodity (+ optional overrides / dateFrom /
/ dateTo). Returns `pnl`validation`meta - validation carries the spread-stationarity gate; meta carries
/ the run record (the R12-shape artifact the workflow audits + attributes).
.template.scs.run:{[inputs]
    cfg:.cfg.strategy.seasonalCalSpread,$[`overrides in key inputs; inputs`overrides; ()!()];
    comm:$[`commodity in key inputs; inputs`commodity; `CRUDE];
    fromD:$[`dateFrom in key inputs; inputs`dateFrom; 0Nd];
    toD:$[`dateTo in key inputs; inputs`dateTo; 0Nd];
    run:.template.scs.__replay[comm;fromD;toD;cfg];
    spread:exec frontReturn from run`steps;                  / not used for stationarity; recompute on the spread
    sp:.template.scs.__spreadSeries[.data.hdb.curveHistory[comm;.data.hdb.dates comm]; cfg`deferredIdx; cfg`rollBuffer];
    stat:.template.rv.stationarity sp`spread;
    pnl:select date:stepDate, pnl:stepPnl from run`steps;
    `pnl`validation`meta!(
        pnl;
        (enlist `stationarity)!enlist `passed`detail!(stat`stationary; "spread-stationarity (seasonal): AR1=",(string stat`ar1),", half-life=",(string stat`halfLife),", stationary=",string stat`stationary);
        `shape`commodity`nDays`ar1`halfLife`edgeSource`runRecord!(
            `seasonalCalSpread; comm; count pnl; stat`ar1; stat`halfLife; `structural; run))
 };

.template.register[`seasonalCalSpread; .template.scs.run;
    `contractVersion`description`shape`composes`in`out!(
        1;
        "seasonally-adjusted crude calendar-spread mean-reversion (the R16 capstone): fade the same-month z, trade the actual legs, emit an auditable run record";
        `seasonalCalSpread;
        `seasonalSpreadReversion;
        (enlist `commodity)!enlist `symbol;
        `pnl`validation`meta!`table`dict`dict)];

-1 "seasonal_calendar_spread.q loaded - .template.scs.* seasonalCalendarSpread template registered (R16 capstone)";
