/ evidence/evidence.q - the Evidence Audit (.evidence.*) (v0.78, Research OS R13)
/ ----------------------------------------------------------------------------
/ ARCHITECTURE.md Part II 11.8(e) / 13-R13. The evidence-layer GATE that closes the loop the door
/ opened: R9's accessor PREVENTS look-ahead by construction; R12's replay RECORDS the provenance
/ (every step's as-of data access, the book, the fills, the rolls, the PnL components); R13 VERIFIES
/ it - a DETERMINISTIC pass/fail audit that bites on a replay run record BEFORE the governance gates
/ ever see its PnL. Wired exactly like carded gating: a HARD precondition, FAIL -> reject, the gates
/ NEVER run, and gov is left UNMODIFIED. After R13 the only PnL the gates judge is PnL PROVEN to
/ correspond to something a desk could actually have traded.
/ -
/ THE AUDIT MUST BITE on every check (an audit that always passes manufactures false confidence): for
/ EVERY check there is a synthetic test that feeds a deliberately-contaminated run and confirms the
/ audit CATCHES it - the same discipline as R2 conformance / R5 .cards.audit / R6 stationarity.
/ -
/ LAYER: ABOVE backtest/ (consumes R12's run record) and BELOW gov/ (reads the zone CONFIG
/ .cfg.gov.zones, NOT gov's logic). Reads state/ (provenance + as-of universe), roll/ (the roll
/ consistency check is on the run record), backtest/ (the run-record shape), config/ (zones + tol) -
/ all downward. Loads after backtest/, before gov/; opens nothing at import. ADDITIVE: it changes no
/ existing compute path, leaves gov UNMODIFIED, and leaves .workflow.run UNCHANGED.
/ -
/ INFRASTRUCTURE (a gate, like the gov gates + the replay engine) - NOT registered as an R2
/ capability and NOT carded (capabilities/templates/strategies get registered + carded; gates do not).
/ -
/ RESERVED-NAME NOTE: parameter `asOf` NEVER `asof`; `comm` NOT `commodity`. Builtins not shadowed
/ (all/any/sum/sums/deltas/differ/min/max/first/last/get/key/value/cols/within).
/ ----------------------------------------------------------------------------

.evidence.defaultConfig:{[] .cfg.evidence};

/ The train+validate boundary date for a commodity, computed from the zone CONFIG .cfg.gov.zones (the
/ SAME floor-cut formula gov uses) + the HDB date axis - a LIGHT cross-check, NOT a call into gov's
/ logic and NOT a second holdout seal (the seal stays in gov, R3b). Returns the last trainValidate date.
.evidence.__trainValidateBoundary:{[comm]
    z:.cfg.gov`zones;
    d:.data.hdb.dates comm;
    n:count d;
    nTV:(floor n*z`trainFrac)+floor n*z`validateFrac;
    d nTV-1
 };

/ ── the deterministic verifier ───────────────────────────────────────────────

/ .evidence.audit[run] -> `pass`checks`reason. Pure function of the run record (+ state/roll/config),
/ deterministic, no RNG. `checks is a dict checkName->boolean; ANY false -> pass=0b with the reason
/ naming the failed checks. `run is an R12 replay run record (`meta`steps`rollEvents`provenance).
.evidence.audit:{[run]
    cfg:.cfg.evidence;
    tol:cfg`reconcileTol;
    steps:run`steps;
    re:run`rollEvents;
    comm:(run`meta)`commodity;
    if[0=count steps; :`pass`checks`reason!(0b;(enlist `nonEmpty)!enlist 0b;"empty run record")];

    / 1. LOOK-AHEAD: every step's asOf is its own date, and no as-of access reached beyond asOf
    / (provDateTo = the max date the step's door slice saw; R9 prevents date>asOf, the audit verifies).
    lookAhead:(all (steps`asOf)=steps`stepDate) and not any (steps`provDateTo)>steps`asOf;

    / 2. UNIVERSE: every FILL's contract was a live (non-expired, expiry>asOf) contract as of its date.
    / (Plain-vector filters, NOT qSQL: select/exec on a LOCAL table inside a lambda throws 'assign.)
    fillMask:0f<>steps`filledQty;
    fillAsOf:(steps`asOf) where fillMask;
    fillContract:(steps`activeContract) where fillMask;
    universe:$[0=count fillAsOf; 1b;
        all {[comm;d;c] sl:(.state.asof[d;comm])`data; expd:first (sl`expiry) where c=sl`contractYM;
             (not null expd) and expd>d}[comm]'[fillAsOf;fillContract]];

    / 3. ROLL-RESPECTED: the recorded rollEvents are exactly the transitions in the active-contract
    / series (no fabricated / missing roll), and every roll is forward in contract (no arbitrary
    / backward roll). The active series itself is R11-produced by the replay; this verifies the record
    / is internally rule-consistent (defense in depth - the roll ENGINE is R11, this checks the trail).
    ac:steps`activeContract;
    chg:(where differ ac) except 0;
    expected:`date xasc ([] date:(steps`stepDate) chg; fromContract:ac chg-1; toContract:ac chg);
    actual:`date xasc ([] date:re`date; fromContract:re`fromContract; toContract:re`toContract);
    rollRespected:(expected~actual) and $[0=count re; 1b; all (re`toContract)>re`fromContract];

    / 4. COSTS-APPLIED: a run with fills carries a positive transaction cost (unless requireCosts off).
    costsApplied:$[cfg`requireCosts; $[0=count fillAsOf; 1b; 0f<sum (steps`proportionalCost) where fillMask]; 1b];

    / 5. FILLS-PRESENT: no position change without a corresponding fill (no position out of nowhere).
    posDelta:deltas steps`position;
    fillsPresent:not any ((posDelta<>0f) and (steps`filledQty)=0f);

    / 6. BOOK-TIES-TO-FILLS: the position book equals the cumulative fills.
    bookTies:tol>max abs (steps`position)-sums steps`filledQty;

    / 7. PnL-TIES: stepPnl == positionPnl (= position x price-move) - totalCost - financingCost, per step.
    resid:(steps`stepPnl)-((steps`positionPnl)-((steps`totalCost)+steps`financingCost));
    pnlTies:tol>max abs resid;

    / 8. DATE-RANGE CONTAINMENT: the run stayed inside train+validate (the holdout was not spanned).
    containment:(max steps`stepDate)<=.evidence.__trainValidateBoundary comm;

    checks:`lookAhead`universe`rollRespected`costsApplied`fillsPresent`bookTies`pnlTies`containment!(
        lookAhead;universe;rollRespected;costsApplied;fillsPresent;bookTies;pnlTies;containment);
    pass:all value checks;
    reason:$[pass; "all evidence checks passed";
        "evidence FAILED: ",", " sv string (key checks) where not value checks];
    `pass`checks`reason!(pass;checks;reason)
 };

/ ── the fail-safe precondition wrapper (mirrors .cards.gatedRun) ──────────────

/ .evidence.gatedRun[hypoId;run;runner;axis] - audit the replay run; on PASS delegate to the gates
/ (.gov.runFull); on FAIL return a reject verdict `evidenceFailed and the gates NEVER run (gov is not
/ invoked). gov is UNMODIFIED; the existing .workflow.run is UNCHANGED. A run that cannot be PROVEN
/ clean is rejected, not judged - mirroring R3b's fail-safe verdicts + R5's carded gating.
.evidence.gatedRun:{[hypoId;run;runner;axis]
    aud:.evidence.audit run;
    if[not aud`pass;
        :enlist `bucket`nObs`verdict`failedGate`tradeable`netSharpe`dsr`postHoc`holdoutPassed`reason!(
            `all;0;`evidenceFailed;`evidence;0b;0n;0n;0b;0b;("REJECTED before gating: ",aud`reason))];
    .gov.runFull[hypoId;runner;axis]
 };

-1 "evidence.q loaded - .evidence.* deterministic evidence audit + fail-safe gatedRun ready (gov unmodified)";
