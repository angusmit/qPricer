/ state/state.q - point-in-time market state + as-of discipline (.state.*) (v0.74, Research OS R9)
/ ----------------------------------------------------------------------------
/ ARCHITECTURE.md Part II 11.8 / 13-R9. The FIRST milestone of the evidence layer and the
/ chokepoint the rest of it stands on: ONE as-of accessor that is the ONLY door to history,
/ plus a first-class Market State object every later foundation milestone (curve R10, roll
/ R11, replay R12, evidence audit R13, season/carry R14, the first real strategy R16) consumes.
/ -
/ WHY ONE DOOR: point-in-time discipline fails when it is a rule every future module must
/ remember. Make it ONE function to audit instead. The accessor PREVENTS look-ahead by
/ construction (it returns only rows knowable as of `asOf`) and logs PROVENANCE; R13's evidence
/ audit later VERIFIES by reading that provenance. Get the door right and the evidence layer
/ inherits point-in-time safety for free.
/ -
/ LOW layer: loads AFTER data/ (it WRAPS the unchanged .data.hdb.* query layer) and BEFORE the
/ consumers (regime/ etc.). Reads data/ downward only; opens no HDB at import (the accessor
/ requires the HDB already opened, like the query layer). ADDITIVE: it does NOT rewire the
/ existing direct readers (regime/, backtest/ keep their current HDB access - rebasing them is
/ a later deferral). Registered as an R2 `state capability; carded via R5.
/ -
/ RESERVED-NAME NOTE: the parameter is `asOf` - NEVER `asof` (the q as-of-join keyword);
/ shadowing it locally is a silent-failure trap.
/ ----------------------------------------------------------------------------

.state.defaultConfig:{[] .cfg.state};

/ ── the provenance log (append-only; what each accessor call returned) ───────

.state.__emptyProv:{[]
    ([] commodity:`symbol$(); asOf:`date$(); dateFrom:`date$(); dateTo:`date$();
        nRows:`long$(); contracts:(); builtAt:`timestamp$())
 };
.state.__initProv:{[] if[not `provLog in key `.state; .state.provLog:.state.__emptyProv[]];};

/ Provenance stamp for an as-of slice (plain dict; no global surprises beyond the append-only log).
.state.__provenance:{[asOf;comm;slice]
    `commodity`asOf`dateFrom`dateTo`nRows`contracts`builtAt!(
        comm; asOf;
        $[count slice; min slice`date; 0Nd];
        $[count slice; max slice`date; 0Nd];
        count slice;
        $[`contractYM in cols slice; asc distinct slice`contractYM; `long$()];
        .z.p)
 };
.state.__logProvenance:{[stamp] .state.__initProv[]; .state.provLog:.state.provLog upsert stamp; stamp};

/ ── the as-of accessor: the ONLY door to history for foundation code ─────────

/ Returns ONLY data knowable as of `asOf` for `comm`, via the existing .data.hdb.* layer
/ (the global `futures` table the query layer loaded). Revision-AWARE by design: if the source
/ carries a `validDate column it takes the latest revision with validDate<=asOf (per observation
/ date); futures carry no revision columns, so they take the simple date<=asOf path (we do NOT
/ fabricate revisable data we do not have). Returns `data`provenance and appends the provenance
/ stamp to the append-only .state.provLog (advisory until R13 consumes it). Parameter `asOf`,
/ never `asof`. (`comm` is the parameter, NOT the column, so the qSQL where-clause is unambiguous.)
.state.asof:{[asOf;comm]
    .data.hdb.__requireLoaded[];
    hasRev:`validDate in cols futures;
    slice:$[hasRev;
        / revisable: knowable observation dates, latest published revision as of asOf
        [pre:select from futures where commodity=comm, date<=asOf, validDate<=asOf;
         select from pre where validDate=(max;validDate) fby date];
        / futures: the simple point-in-time filter
        select from futures where commodity=comm, date<=asOf];
    stamp:.state.__provenance[asOf;comm;slice];
    .state.__logProvenance stamp;
    `data`provenance!(slice;stamp)
 };

/ ── the Market State object: the contract signals/templates consume ──────────

/ Build the point-in-time Market State for (asOf, comm), assembled lazily from the accessor
/ (built per call, NOT materialised for all dates). A flat dict of typed fields with nested
/ sub-objects: asOf, commodity, effectiveDate (the latest trading date <= asOf), curve (the
/ as-of curve, the SAME shape .data.hdb.curveAt returns so it is consistent with what regime/
/ sees), spreads (adjacent calendar spreads from the curve), features (a placeholder the curve
/ engine R10 + season/carry R14 fill - declared now so the contract is stable), universe (the
/ tradable contracts NOT expired as of asOf), refs (cost model + roll-rule handles; the real
/ roll map is R11), and provenance. Pure given (asOf, comm) + the HDB.
.state.build:{[asOf;comm]
    a:.state.asof[asOf;comm];
    slice:a`data;
    if[0=count slice; '"state.build: no data for ",(string comm)," as of ",string asOf];
    effDate:max slice`date;
    curve:.data.hdb.curveAt[comm;effDate];
    spreads:$[1<count curve;
        ([] nearTenor:(-1_curve`tenor); farTenor:(1_curve`tenor);
            nearContract:(-1_curve`contractYM); farContract:(1_curve`contractYM);
            spread:(-1_curve`price)-1_curve`price);
        ([] nearTenor:`float$(); farTenor:`float$(); nearContract:`long$(); farContract:`long$(); spread:`float$())];
    cfg:.cfg.state;
    liveExpiry:$[cfg`expiryStrictlyAfter; asOf<; asOf<=];
    universe:asc distinct exec contractYM from slice where liveExpiry expiry;
    refs:`costModel`rollRule`commodity!(`dailyFillCost; comm; comm);
    `asOf`commodity`effectiveDate`curve`spreads`features`universe`refs`provenance!(
        asOf; comm; effDate; curve; spreads; (`symbol$())!(); universe; refs; a`provenance)
 };

/ ── the no-look-ahead invariants (tests assert them; R13 will reuse them) ────

/ (i) the accessor never returns a row with date>asOf.
.state.invariant.asofRespected:{[asOf;data] not any data[`date]>asOf};
/ (ii) the as-of live curve carries no expired contract (every curve expiry is on/after the
/ effective date) - so the tradable universe built from it excludes contracts expired by asOf.
.state.invariant.universeLive:{[st] all (st`effectiveDate) <= (st`curve)`expiry};
/ (iv) the provenance stamp matches exactly what was returned.
.state.invariant.provenanceConsistent:{[a]
    d:a`data; p:a`provenance;
    ((p`nRows)=count d) and ((p`dateTo)~$[count d; max d`date; 0Nd]) and (p`dateFrom)~$[count d; min d`date; 0Nd]};

/ ── register as an R2 `state capability ──────────────────────────────────────

.state.contract:`version`requiredIn`requiredOut!(
    1;
    `asOf`commodity!`date`symbol;
    `asOf`commodity`curve`universe`provenance!`date`symbol`table`symbol`dict);
.registry.new[`state; .state.contract];
.state.register:{[name;fn;manifest] .registry.register[`state;name;fn;manifest]};
.state.get:{[name] .registry.get[`state;name]};
.state.list:{[] .registry.list `state};
.state.conforms:{[name] .registry.conforms[`state;name]};

.state.register[`marketState; .state.build;
    `contractVersion`description`in`out!(
        1;
        "point-in-time Market State for (asOf, commodity): the single door to history + the state object";
        `asOf`commodity!`date`symbol;
        `asOf`commodity`curve`universe`provenance!`date`symbol`table`symbol`dict)];

.state.__initProv[];
-1 "state.q loaded - .state.* as-of accessor + Market State ready (the single door to history; not opened)";
