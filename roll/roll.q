/ roll/roll.q - roll discipline (.roll.*) (v0.76, Research OS R11)
/ ----------------------------------------------------------------------------
/ ARCHITECTURE.md Part II 11.8(c) / 13-R11. Commodity backtests are most often wrong because
/ of roll mistakes - rolling on data the strategy could not have seen, or backtesting an
/ abstract continuous series instead of the contracts you would actually trade. R11 makes roll
/ discipline a config + a single function: which contract is ACTIVE as of a date, decided
/ DETERMINISTICALLY FROM AS-OF DATA ONLY (through R9's door), with roll events recorded and the
/ continuous series demoted to an analytics-only derived view.
/ -
/ THE PRINCIPLE: trade actual contracts (the roll map names them); never trade or signal off the
/ continuous series. LOW layer: depends ONLY on R9's door (.state.asof) for as-of expiry/volume/
/ prices; independent of curve/ (placed after it for tidiness). Reads state/ (+ data/) downward;
/ opens nothing at import. ADDITIVE: it does NOT rewire the existing readers.
/ -
/ HONEST DATA SCOPE: the HDB carries VOLUME + EXPIRY but NO open interest, so oi_switch is
/ RECOGNISED but FLAGGED unavailable (we do not fabricate OI) - the same way R9 designed
/ revision-awareness without fabricating revision columns.
/ -
/ RESERVED-NAME NOTE: parameter `asOf` NEVER `asof`; `comm` NOT `commodity` (the qSQL self-compare
/ trap). Builtins not shadowed (next/prev/ratios/deltas/differ/first/last/bin/binr).
/ ----------------------------------------------------------------------------

.roll.defaultConfig:{[] .cfg.rolls};

/ The per-commodity rule (falls back to the `default rule).
.roll.__rule:{[comm] $[comm in key .cfg.rolls; .cfg.rolls comm; .cfg.rolls`default]};

/ The live contracts as of asOf (expiry>asOf), with their trailing-window volume, sorted by
/ expiry - obtained THROUGH R9's door (date<=asOf), so the roll decision sees no future data.
.roll.__liveAsOf:{[asOf;comm;volWin]
    slice:(.state.asof[asOf;comm])`data;
    if[0=count slice; '"roll: no as-of data for ",string comm];
    effDate:max slice`date;
    `expiry xasc 0!select expiry:first expiry, vol:sum volume, effDate:max date
        by contractYM from select from slice where expiry>asOf, date within (effDate-volWin;effDate)
 };

/ ── rule dispatchers (each returns `active`weights`reason; as-of-only) ────────

/ days_before_expiry (default): hold the nearest-expiry contract whose (expiry-asOf) > rollDays
/ (you have rolled out of any nearer one). In the W days before that contract's own roll point
/ (heldD2e in (rollDays, rollDays+W]) blend held->next, so the roll is not a single-day jump.
.roll.__daysBeforeExpiry:{[asOf;live;cfg]
    rd:cfg`rollDays; w:cfg`W;
    d2e:`int$(live`expiry)-asOf;
    idxs:where d2e>rd;
    heldIdx:$[count idxs; first idxs; (count live)-1];
    held:(live`contractYM) heldIdx;
    heldD2e:d2e heldIdx;
    hasNext:heldIdx<(count live)-1;
    nextC:$[hasNext; (live`contractYM) heldIdx+1; held];
    $[hasNext and heldD2e<=rd+w;
        [frac:(`float$heldD2e-rd)%w; `active`weights`reason!($[frac>=0.5;held;nextC]; (held;nextC)!(frac;1f-frac); `rolling)];
        `active`weights`reason!(held; (enlist held)!enlist 1f; `held)]
 };

/ volume_switch: front vs next by trailing as-of volume; roll when next's volume exceeds front's.
.roll.__volumeSwitch:{[asOf;live;cfg]
    front:(live`contractYM) 0;
    $[1<count live;
        [nextC:(live`contractYM) 1; rollIt:(live[`vol] 1)>live[`vol] 0;
         `active`weights`reason!($[rollIt;nextC;front]; (enlist $[rollIt;nextC;front])!enlist 1f; $[rollIt;`volumeRolled;`held])];
        `active`weights`reason!(front; (enlist front)!enlist 1f; `held)]
 };

/ fixed_calendar: roll to next on/after the rollDays-th day of the front's expiry month.
.roll.__fixedCalendar:{[asOf;live;cfg]
    front:(live`contractYM) 0;
    rollDate:(`date$`month$(live`expiry) 0)+cfg`rollDays;
    $[(1<count live) and asOf>=rollDate;
        [nextC:(live`contractYM) 1; `active`weights`reason!(nextC; (enlist nextC)!enlist 1f; `calendarRolled)];
        `active`weights`reason!(front; (enlist front)!enlist 1f; `held)]
 };

/ ── the roll-rule engine (.roll.active) ──────────────────────────────────────

/ The active contract(s) as of asOf, decided from the door's as-of data ONLY. Returns
/ `active (the primary contractYM), `weights (contractYM!weight - one contract normally, two
/ during a roll window), `reason, plus asOf/commodity/effectiveDate/ruleType.
.roll.active:{[asOf;comm]
    cfg:.roll.__rule comm;
    if[(cfg`type)=`oi_switch;
        '"roll.active: oi_switch is recognised but UNAVAILABLE (the HDB has no open-interest column)"];
    live:.roll.__liveAsOf[asOf;comm;cfg`W];
    if[0=count live; '"roll.active: no live contract for ",(string comm)," as of ",string asOf];
    res:$[(cfg`type)=`days_before_expiry; .roll.__daysBeforeExpiry[asOf;live;cfg];
          (cfg`type)=`volume_switch; .roll.__volumeSwitch[asOf;live;cfg];
          (cfg`type)=`fixed_calendar; .roll.__fixedCalendar[asOf;live;cfg];
          '"roll.active: unknown rule type ",string cfg`type];
    res,`asOf`commodity`effectiveDate`ruleType!(asOf;comm;max live`effDate;cfg`type)
 };

/ ── roll events (splayed, UNPARTITIONED rollEvents) ──────────────────────────

.roll.__emptyEvents:{[]
    ([] date:`date$(); commodity:`symbol$(); fromContract:`long$(); toContract:`long$(); reason:`symbol$())
 };

/ The active contract+reason as of a date, or the generic null on a date with no live contract
/ (the data edge: the last bars are often expiry-day artifacts with no un-expired deferred
/ contract). Series builders below drop those dates - the series is defined only where it trades.
.roll.__activeOrNull:{[comm;d] @[{[comm;d] r:.roll.active[d;comm]; `active`reason!(r`active;r`reason)}[comm;];d;{[e] (::)}]};

/ Walk the dates and emit a roll-event row whenever the active contract changes.
.roll.events:{[comm;fromDate;toDate]
    dts:.data.hdb.dates comm;
    dts:dts where dts within (fromDate;toDate);
    if[0=count dts; :.roll.__emptyEvents[]];
    recs:.roll.__activeOrNull[comm] each dts;
    live:not (::)~/:recs;                 / drop dates with no live contract (data-edge artifacts)
    dts:dts where live; recs:recs where live;
    if[0=count dts; :.roll.__emptyEvents[]];
    actives:recs`active; reasons:recs`reason;
    prevA:prev actives;
    chg:where (not null prevA) and actives<>prevA;
    ([] date:dts chg; commodity:comm; fromContract:prevA chg; toContract:actives chg; reason:reasons chg)
 };

/ Persist / load the rollEvents splay (splayed, UNPARTITIONED, p# on commodity, .Q.en + set -
/ the SAME storage decision as the HDB; NOT date-partitioned). gitignored; tests synthetic only.
.roll.eventsPersist:{[hdbPath;events]
    t:.Q.en[hsym `$hdbPath; update `p#commodity from `commodity`date xasc events];
    (hsym `$hdbPath,"/rollEvents/") set t;
    `rows!enlist count t
 };
.roll.eventsOpen:{[hdbPath]
    sp:hsym `$hdbPath,"/rollEvents/";
    if[0=count key sp; '"roll.eventsOpen: no rollEvents at ",hdbPath];
    get sp
 };

/ ── the continuous-adjusted view (ANALYTICS-ONLY, NON-POINT-IN-TIME) ─────────

/ A DERIVED back-adjusted continuous price series: the active contract's settle per date,
/ spliced at each roll using the new contract's OWN move (difference method) so the series has no
/ gap. WARNING: back-adjustment rewrites historical levels at every roll, so this series is
/ NON-POINT-IN-TIME and is for CHARTING / long-run vol ONLY - never trade or signal off it (the
/ point-in-time-safe, tradable path is .roll.active, the actual contracts). R13 flags any strategy
/ that signals off this series. (The active-contract CHOICE per date is still as-of-clean via
/ .roll.active; only the price back-adjustment is retrospective.)
.roll.continuous:{[comm]
    dts:.data.hdb.dates comm;
    ftk:`date`contractYM xkey select date,contractYM,settle from futures where commodity=comm;
    recs:.roll.__activeOrNull[comm] each dts;
    keep:not (::)~/:recs;                 / defined only where a live active contract exists
    dts:dts where keep;
    if[0=count dts; '"roll.continuous: no live active contract on any date for ",string comm];
    active:(recs where keep)`active;
    lkp:{[ftk;ds;cs] (ftk ([] date:ds; contractYM:cs))`settle};
    raw:lkp[ftk;dts;active];
    rolled:active<>prev active; rolled[0]:0b;
    / within-contract day return: on a roll, use the NEW contract's own move (its prior-day settle).
    newPrior:lkp[ftk; prev dts; active];
    ret:?[rolled; raw-newPrior; raw-prev raw];
    ret:0f^ret; ret[0]:0f;
    ([] date:dts; contract:active; rawFront:raw; adjusted:(first raw)+sums ret; rolled)
 };

/ ── register as an R2 `roll capability ───────────────────────────────────────

.roll.contract:`version`requiredIn`requiredOut!(
    1;
    `asOf`commodity!`date`symbol;
    `active`weights`reason!`long`dict`symbol);
.registry.new[`roll; .roll.contract];
.roll.register:{[name;fn;manifest] .registry.register[`roll;name;fn;manifest]};
.roll.get:{[name] .registry.get[`roll;name]};
.roll.list:{[] .registry.list `roll};
.roll.conforms:{[name] .registry.conforms[`roll;name]};

.roll.register[`rollEngine; .roll.active;
    `contractVersion`description`in`out!(
        1;
        "as-of-only active-contract roll rules (days_before_expiry/volume_switch/fixed_calendar; oi_switch flagged unavailable) + window blend";
        `asOf`commodity!`date`symbol;
        `active`weights`reason!`long`dict`symbol)];

-1 "roll.q loaded - .roll.* roll discipline ready (as-of-only active contract + events + continuous view)";
