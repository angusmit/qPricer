# cards/ — Model cards: the knowledge plug-in (ARCHITECTURE.md Part II knowledge plug-ins, R5)

## Purpose
A model card is a structured, queryable record per capability that **synthesises** the rest of the system: what the capability **is** (the R2 contract it satisfies), what it **assumes**, which of the three **edge sources** it claims (for signals/strategies), which **regimes** it's valid in, its linked **risk memory** (R4), and — the part that makes a card honest rather than marketing — its **validation status derived from the gov ledger** (R3/R3b), never asserted. It turns "trust me, this works" into a record you can query and that can't overstate itself.

## Dependencies
A **HIGH** layer — a card sits *above* the things it documents. It reads (all legal downward): the R2 registries (`.model`/`.signal`/`.calibrator`/`.execution.fillModel.list` + `.strategy.registeredStrategies`), the gov ledger + hypotheses (`.gov.trials` / `.gov.hypoTbl`, recomputing the deflated Sharpe via `.gov.deflatedSharpe`), and the R4 regime library (`.regime.library`). It **must not** be imported by `gov/`/`regime/`/`backtest/`. Loads after `gov/`, before `apps/`. **Never opens the HDB at import** (`.cards.open` uses `get`, not `\l`).

## Modules
- `cards.q` — `.cards.*`: the card schema/store, lookup, the gov-derived validation status, the regime risk-memory link, and the audit.

## Key API
- `.cards.list[]` / `.cards.get[capName]` / `.cards.forCapability[kind;name]` — card lookup. (Parameters are `capName`/`capKind`, deliberately **not** the column names `capabilityName`/`capabilityKind` — a qSQL `where capabilityName=capabilityName` compares the column to itself.)
- `.cards.validationStatus[capName]` → a record `status`/`nTrials`/`headlineNetSharpe`/`headlineDsr`/`deflationPass`/`holdoutUsed`/`tradeable`, **derived** from the gov ledger via the card's `govHypoId`. `ungated` when nothing is logged; `inResearch` (logged, no holdout); `validated`/`rejected` (holdout passed/failed). `tradeable` iff the sealed holdout passed — it can never read "validated" when the ledger says otherwise.
- `.cards.riskMemory[capName]` → the linked regime-library failure modes (reuses R4, via the card's `riskMemoryKey`).
- `.cards.audit[]` → the honesty check (below). `.cards.buildTable[hdbPath]` / `.cards.open[hdbPath]` — persist/load the `modelCards` splay (via `get`/`set` + `.Q.en`).

## The schema
`modelCards` (curated in `.cfg.cards`, persisted by `scripts/build_model_cards.q`): `cardId`, `capabilityKind`, `capabilityName` (matches the registry name), `version`, `intendedUse`, `assumptions`, `edgeSource` (`riskPremium`/`structural`/`informational` for signals/strategies, `na` for pricers/fill models), `regimeApplicability`, `riskMemoryKey` (→ the regime library, `na` if none), `govHypoId` (→ the gov record, `` ` `` if not gated), `owner`, `asOf`. The full narrative per card is in `docs/MODEL_CARDS.md`.

## The audit (the honesty check that CAN FAIL)
`.cards.audit[]` (pure core `.cards.__auditAgainst`) walks **every** registered capability (the four R2 kinds + the 24 strategies) and flags: an **undocumented** registered capability; a card missing a required section (`intendedUse`/`assumptions`/`regimeApplicability` empty); a **signal/strategy card with no named edge source** ("name which of the three or it's noise"); and an **orphan card** (capabilityName not registered). It passes only fully-documented, edge-named, resolvable cards — the same discipline as the gov gates and the R2 conformance check. Today 7 capabilities are carded; the audit honestly reports the rest as a coverage gap (curate, don't auto-spam).

## Notes
- **Validation is derived, never asserted** — the crux of an honest card. The status reflects the actual governance record; promote a capability by passing the gates, not by editing a card.
- ADDITIVE: reads existing artifacts, changes no compute path, so the full suite stays byte-identical. Demo: `apps/examples/model_card_show.q` (one screen: contract + assumptions + edge + live derived validation + risk memory + audit). Build: `scripts/build_model_cards.q`.
- The knowledge plug-in R6 (problem-templates) and R7 (agents) read from when deciding what to run and how much to trust it.
