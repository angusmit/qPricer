# Capability contracts — the plug-in spine

**Version:** 0.68 (Research OS R2) · **Mechanism:** `core/registry.q` (`.registry.*` + the per-kind
registries) · **Populated by:** `core/registry_populate.q` · **Inventory:** `apps/examples/registry_inventory.q`.

The whole architecture rests on **one small generic spine + plug-ins through fixed contracts**. A new
capability does not edit the engine; it **registers** itself against a versioned contract and is admitted
only if it **conforms**. The conformance check **can fail** — a plug-in that does not honour its contract
is rejected, not rubber-stamped (the same discipline as the gov gates and the regime library: a check that
can't say no is theatre).

## The registry mechanism (`.registry.*`)
- `.registry.new[kind;contract]` — declare a capability kind + its contract (idempotent).
- `.registry.register[kind;name;fn;manifest]` — register a plug-in: `fn` is a **handle to the unchanged
  capability function**; `manifest` is its declared I/O schema. Unknown kind → signals; **duplicate name →
  overwrite-with-warn** (latest wins).
- `.registry.get[kind;name]` → `` `name`kind`fn`manifest `` (unknown → signals). `.registry.list[kind]` →
  the registered names. `.registry.unregister[kind;name]` / `.registry.dropKind[kind]`.
- `.registry.conforms[kind;name]` → bool. `.contracts.verify[]` → a table (kind, name, pass, reason) over
  **every** registered capability — the spine-wide health check.

Per-kind registries are thin specialisations carrying their contract: `.model.*`, `.signal.*`,
`.calibrator.*`, `.execution.fillModel.*` (each with `register`/`get`/`list`/`conforms`).

## Manifests and conformance
A **manifest** is a dict ``  `contractVersion`in`out  `` (+ optional `description`), where `in`/`out` are
**schemas** — dicts of `fieldName → typeTag` (`` `float`long`symbol`boolean`char`dict`table ``). A plug-in
**conforms** to its kind's contract iff:
1. the manifest's `contractVersion` matches the contract's `version`, AND
2. every contract-required **input** key is present in the manifest `in` with a matching type tag, AND
3. every contract-required **output** key is present in the manifest `out` with a matching type tag.

The manifest **may declare extra fields** — a superset conforms. A manifest missing a required key, or
declaring a wrong type tag, or targeting an incompatible contract version, **fails** (proven in
`tests/registry/test_contracts_conform.q`).

## The contracts (v1) and the registered capabilities

| Kind | Contract — required in | Contract — required out | Registered (handle) |
|---|---|---|---|
| `model` (pricer) | `spot`,`strike`,`expiry`,`rate`,`vol` (float) | `unitPrice`,`notionalPrice` (float) | `blackScholesFdm` → `.engine.priceOption` |
| `signal` | `curveHistory` (table), `sigCfg` (dict) | `path` (table) | `commoditySignalPath` → `.strategy.path.commoditySignals` |
| `calibrator` | `marketCurve` (table), `modelType` (symbol), `calCfg` (dict) | `calibratedParams` (dict), `fitRmse` (float), `status` (symbol) | `schwartz2Curve` → `.commodity.calibrateCurve` |
| `fillModel` (execution) | `order` (float), `ctx` (dict), `cfg` (dict) | `fillPrice` (float), `totalCost` (float) | `dailyFillCost` → `.exec.fill` |

All four registered capabilities pass `.contracts.verify[]`.

## Scope + notes
- **Bounded (R2):** the mechanism + the contracts + registering what exists + a conformance test. It is
  NOT a migration of every call site and NOT a new capability. **Existing direct call paths are untouched**
  (the backtest still finds a strategy and the pricer is still invoked exactly as before); the registry is
  the canonical lookup for **new** work and for introspection/governance. Registration is **metadata only**
  — no compute path changed, so the suite stays byte-identical.
- **Strategies** keep their own existing registry (`.strategy.register` / `.strategy.registry`),
  deliberately **not rebased** (preserving byte-identity); the four contracts above are the *new* kinds.
- This is the foundation **R5** (model cards — a card registry the gov Gate 0 reads), **R6** (problem
  templates), and **R7** (agents) plug into: each declares a contract and snaps in.
- **Adding a plug-in:** call `.<kind>.register[name; fn; manifest]` with a manifest that declares the
  contract's required keys/types (plus whatever extra it produces), then confirm `.<kind>.conforms[name]`.
