/ carry/carry.q - carry / storage economics (.carry.*) (v0.79, Research OS R14)
/ ----------------------------------------------------------------------------
/ ARCHITECTURE.md Part II 11.8(f) / 13-R14. An evidence-layer FEATURE capability: the CARRY / storage
/ economics the carry strategies are built on - implied carry, convenience yield, cash-and-carry fair
/ value, the carry signal, an inventory-tightness proxy. Derived from R10's curve (through R9's door),
/ so it is point-in-time by construction.
/ -
/ HONEST DATA SCOPE: the implied carry comes from the CURVE (we have it - it reuses R10's annualised
/ rollYield). The convenience yield + cash-and-carry fair value need a rate r + a storage cost, which we
/ SUPPLY from .cfg.carry (config, NOT market-observed) - so those outputs are ASSUMPTION-DEPENDENT and
/ flagged (`assumptionDependent). The inventory-tightness is a PROXY from the degree of backwardation;
/ we have NO real inventory data, so it is labelled a proxy and inventory is never fabricated (the same
/ discipline as R11's no-open-interest).
/ -
/ SIGN CONVENTION (matching the milestone spec + the test): impliedCarry is the curve-implied carry/roll,
/ annualised ln(Ffront/Fdeferred)/T - POSITIVE in backwardation, NEGATIVE in contango (reuses R10's
/ rollYield exactly). convenienceYield = (r+storage) - impliedCarry; carrySignal = impliedCarry -
/ (r+storage) = the roll yield net of the assumed cost of carry (the edge signal, flips with the curve).
/ -
/ LOW evidence-tier: loads AFTER curve/ (reuses .curve.build) and reads state/ via it, downward; opens
/ nothing at import. ADDITIVE: reads + computes, edits nothing. Registered as an R2 `carry capability; carded.
/ -
/ RESERVED-NAME NOTE: `asOf` NEVER `asof`; `comm` NOT `commodity`; builtins not shadowed (exp/log/sum/
/ avg/first/last). Mind right-associativity: front*exp costOfCarry*dt = front*(exp(costOfCarry*dt)).
/ ----------------------------------------------------------------------------

.carry.defaultConfig:{[] .cfg.carry};

/ The carry feature dict for (asOf, comm), from R10's curve + the .cfg.carry assumptions.
.carry.features:{[asOf;comm]
    cfg:.cfg.carry;
    r:cfg`riskFreeRate; storage:cfg`storageCost;
    costOfCarry:r+storage;
    cv:.curve.build[asOf;comm];
    curve:cv`curve; feat:cv`features;
    prices:`float$curve`price; tenors:`float$curve`tenor;
    n:count prices;
    front:first prices;
    defIdx:(.cfg.curve`deferredIdx)&n-1;            / the SAME deferred leg the curve's rollYield used
    deferred:prices defIdx;
    dt:(tenors defIdx)-first tenors;                / the front->deferred tenor gap (years)
    rollYield:feat`rollYield;                       / annualised ln(Ffront/Fdeferred)/dt (R10; +ve backwardation)
    impliedCarry:rollYield;                         / the curve-implied carry/roll (+ve backwardation)
    `asOf`commodity`effectiveDate`front`deferred`tenorGap`classification`impliedCarry`costOfCarry`convenienceYield`fairValueCarry`actualDeferred`carrySignal`inventoryTightness`assumptionDependent!(
        asOf; comm; cv`effectiveDate; front; deferred; dt; feat`classification;
        impliedCarry;
        costOfCarry;
        $[null impliedCarry; 0n; costOfCarry-impliedCarry];                 / convenience yield = r+storage - carry
        $[(front>0f) and dt>0f; front*exp costOfCarry*dt; 0n];             / cash-and-carry fair-value deferred F
        deferred;                                                          / the ACTUAL deferred F (compare to fair value)
        $[null impliedCarry; 0n; impliedCarry-costOfCarry];                / carry signal = roll yield net of carry cost
        impliedCarry;                                                      / inventory-tightness PROXY (backwardation degree)
        `convenienceYield`fairValueCarry)                                  / the assumption-dependent outputs
 };

/ ── register as an R2 `carry capability ──────────────────────────────────────

.carry.contract:`version`requiredIn`requiredOut!(
    1;
    `asOf`commodity!`date`symbol;
    `asOf`commodity`impliedCarry`convenienceYield`carrySignal!`date`symbol`float`float`float);
.registry.new[`carry; .carry.contract];
.carry.register:{[name;fn;manifest] .registry.register[`carry;name;fn;manifest]};
.carry.get:{[name] .registry.get[`carry;name]};
.carry.list:{[] .registry.list `carry};
.carry.conforms:{[name] .registry.conforms[`carry;name]};

.carry.register[`carryEconomics; .carry.features;
    `contractVersion`description`in`out!(
        1;
        "carry/storage economics from the curve: implied carry, convenience yield + cash-and-carry fair value (assumption-dependent on r+storage), carry signal, inventory-tightness proxy";
        `asOf`commodity!`date`symbol;
        `asOf`commodity`impliedCarry`convenienceYield`carrySignal!`date`symbol`float`float`float)];

-1 "carry.q loaded - .carry.* carry/storage economics ready (convenience yield + fair value assumption-dependent; inventory a proxy; not opened)";
