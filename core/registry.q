/ core/registry.q - the plug-in spine: a generic registry + per-kind registries +
/ contracts + a conformance check (.registry.* / .model|.signal|.calibrator|
/ .execution.fillModel.* / .contracts.*) (v0.68, Research OS R2)
/ ----------------------------------------------------------------------------
/ ARCHITECTURE.md Part II (integration principle) / 13-R2. The architecture rests on
/ ONE small generic spine + plug-ins through fixed CONTRACTS. This module is that spine:
/   1. a GENERIC registry (.registry.*): new / register / get / list / conforms,
/   2. PER-KIND registries built on it (.model / .signal / .calibrator /
/      .execution.fillModel), each carrying its CONTRACT (declared in/out schema + version),
/   3. a conformance verifier (.registry.conforms / .contracts.verify) that CAN FAIL -
/      a non-conforming plug-in is rejected, not rubber-stamped (same discipline as the
/      gov gates: a check that can't fail is theatre).
/ -
/ LOW module: loads right after config/, before models/ - it depends on nothing above
/ base q. Registration is METADATA ONLY: a plug-in is registered with a manifest (its
/ declared I/O schema) + a HANDLE to the UNCHANGED function. The registry NEVER changes
/ how a capability computes; existing direct call paths are untouched. It is the canonical
/ lookup for new work + introspection/governance, and the foundation R5/R6/R7 plug into.
/ -
/ A manifest is a dict: `contractVersion`in`out (+ optional `description), where `in/`out
/ are schemas - dicts of fieldName -> typeTag (`float`long`symbol`boolean`char`dict`table).
/ Conformance checks the manifest against the kind's contract: version match + every
/ contract-required in/out key present in the manifest with a matching type tag (the
/ manifest MAY declare more - a superset conforms). A deficient manifest FAILS.
/ ----------------------------------------------------------------------------

/ ── generic registry state + primitives ──────────────────────────────────────

.registry.__contracts:()!();   / kind (sym) -> contract dict
.registry.__store:()!();        / kind (sym) -> (name (sym) -> entry dict)

.registry.__assertKind:{[kind]
    if[not kind in key .registry.__contracts; '"registry: unknown kind ",string kind];
 };

/ Create a registry for a capability kind, carrying its contract. Idempotent: re-declaring
/ a kind updates its contract but preserves any already-registered entries.
.registry.new:{[kind;contract]
    .registry.__contracts[kind]:contract;
    if[not kind in key .registry.__store; .registry.__store[kind]:()!()];
    kind
 };

.registry.kinds:{[] key .registry.__contracts};
.registry.contract:{[kind] .registry.__assertKind kind; .registry.__contracts kind};

/ Register a plug-in: fn is a HANDLE to the (unchanged) capability function; manifest is
/ its declared I/O schema (+ metadata). Unknown kind -> signal. Duplicate name ->
/ overwrite-with-warn (the latest registration wins; a warning is emitted).
.registry.register:{[kind;name;fn;manifest]
    .registry.__assertKind kind;
    sub:.registry.__store kind;
    if[name in key sub; -1 "registry.register: overwriting ",(string kind),"/",string name];
    entry:`name`kind`fn`manifest!(name;kind;fn;manifest);
    .registry.__store[kind]:sub,(enlist name)!enlist entry;
    name
 };

/ Remove a single plug-in / a whole kind (used by tests to keep the global registry clean).
.registry.unregister:{[kind;name] .registry.__assertKind kind; .registry.__store[kind]:(.registry.__store kind) _ name; name};
.registry.dropKind:{[kind] .registry.__contracts:.registry.__contracts _ kind; .registry.__store:.registry.__store _ kind; kind};

/ Get a registered entry (`name`kind`fn`manifest). Unknown kind / unknown name -> signal.
.registry.get:{[kind;name]
    .registry.__assertKind kind;
    sub:.registry.__store kind;
    if[not name in key sub; '"registry.get: ",(string kind),"/",(string name)," not registered"];
    sub name
 };

/ The registered names for a kind.
.registry.list:{[kind] .registry.__assertKind kind; key .registry.__store kind};

/ ── conformance (the check that CAN FAIL) ────────────────────────────────────

/ Check a declared schema against a required schema (both dicts fieldName -> typeTag).
/ Returns (ok; missingKeys; badTypeKeys): every required key must be present with a
/ matching type tag (the declared schema MAY carry extra keys - a superset is fine).
.registry.__conformSchema:{[required;declared]
    reqKeys:key required;
    missing:reqKeys where not reqKeys in key declared;
    present:reqKeys where reqKeys in key declared;
    badType:present where not (required present)=declared present;
    (0=count missing,badType; missing; badType)
 };

/ Conformance detail for one registered plug-in: `kind`name`pass`reason.
.registry.__check:{[kind;name]
    .registry.__assertKind kind;
    c:.registry.__contracts kind;
    sub:.registry.__store kind;
    if[not name in key sub; :`kind`name`pass`reason!(kind;name;0b;"not registered")];
    m:(sub name)`manifest;
    mv:$[`contractVersion in key m; m`contractVersion; 0N];
    if[not (c`version)~mv;
        :`kind`name`pass`reason!(kind;name;0b;"contract version mismatch: contract ",(string c`version)," vs manifest ",string mv)];
    inChk:.registry.__conformSchema[c`requiredIn; $[`in in key m; m`in; ()!()]];
    outChk:.registry.__conformSchema[c`requiredOut; $[`out in key m; m`out; ()!()]];
    ok:(inChk 0) and outChk 0;
    reason:$[ok; "conforms to ",(string kind)," contract v",string c`version;
        "non-conforming - in(missing: ",("," sv string inChk 1),"; badType: ",("," sv string inChk 2),") out(missing: ",("," sv string outChk 1),"; badType: ",("," sv string outChk 2),")"];
    `kind`name`pass`reason!(kind;name;ok;reason)
 };

/ Does a registered plug-in conform to its kind's contract? (bool)
.registry.conforms:{[kind;name] (.registry.__check[kind;name])`pass};

/ Run conformance over EVERY registered capability of EVERY kind -> a table
/ (kind, name, pass, reason). The honest spine-wide health check.
.contracts.verify:{[]
    kinds:.registry.kinds[];
    rows:raze {[k] .registry.__check[k;] each .registry.list k} each kinds;
    $[0=count rows;
        ([] kind:`symbol$(); name:`symbol$(); pass:`boolean$(); reason:());
        rows]
 };

/ ── per-kind registries (thin specialisations carrying their contract) ───────

/ model / pricer: in = market state (spot/strike/expiry/rate/vol); out = price.
.model.contract:`version`requiredIn`requiredOut!(
    1;
    `spot`strike`expiry`rate`vol!`float`float`float`float`float;
    `unitPrice`notionalPrice!`float`float);
.registry.new[`model; .model.contract];
.model.register:{[name;fn;manifest] .registry.register[`model;name;fn;manifest]};
.model.get:{[name] .registry.get[`model;name]};
.model.list:{[] .registry.list `model};
.model.conforms:{[name] .registry.conforms[`model;name]};

/ signal: in = a curve history (HDB-derived table) + a signal config; out = a per-date
/ signal/position series (the `path table).
.signal.contract:`version`requiredIn`requiredOut!(
    1;
    `curveHistory`sigCfg!`table`dict;
    (enlist `path)!enlist `table);
.registry.new[`signal; .signal.contract];
.signal.register:{[name;fn;manifest] .registry.register[`signal;name;fn;manifest]};
.signal.get:{[name] .registry.get[`signal;name]};
.signal.list:{[] .registry.list `signal};
.signal.conforms:{[name] .registry.conforms[`signal;name]};

/ calibrator: in = a market curve + a model target + a cal config; out = fitted params +
/ a fit diagnostic + a status.
.calibrator.contract:`version`requiredIn`requiredOut!(
    1;
    `marketCurve`modelType`calCfg!`table`symbol`dict;
    `calibratedParams`fitRmse`status!`dict`float`symbol);
.registry.new[`calibrator; .calibrator.contract];
.calibrator.register:{[name;fn;manifest] .registry.register[`calibrator;name;fn;manifest]};
.calibrator.get:{[name] .registry.get[`calibrator;name]};
.calibrator.list:{[] .registry.list `calibrator};
.calibrator.conforms:{[name] .registry.conforms[`calibrator;name]};

/ fillModel (execution): in = a signed order + market context + a cost config; out = a
/ fill price + a total cost.
.execution.fillModel.contract:`version`requiredIn`requiredOut!(
    1;
    `order`ctx`cfg!`float`dict`dict;
    `fillPrice`totalCost!`float`float);
.registry.new[`fillModel; .execution.fillModel.contract];
.execution.fillModel.register:{[name;fn;manifest] .registry.register[`fillModel;name;fn;manifest]};
.execution.fillModel.get:{[name] .registry.get[`fillModel;name]};
.execution.fillModel.list:{[] .registry.list `fillModel};
.execution.fillModel.conforms:{[name] .registry.conforms[`fillModel;name]};

-1 "registry.q loaded - .registry.* spine + model/signal/calibrator/fillModel contracts ready";
