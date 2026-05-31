# Regime Library — named historical episodes + risk memory

**Version:** 0.67 (Research OS R4) · **Backed by:** the `regimeEpisodes` table (built by
`scripts/build_regime_library.q` from the `regimes` table) + the analogue engine `.regime.analogue.*`.

This is the project's **digested history**: a curated, versioned set of named crude regime episodes,
each with its drivers and — the part that matters — its **risk memory** (the known strategy failure
modes, "what killed people here"), written down so the workflow can ask *"this state resembles what, and
what killed strategies in that state?"* from an explicit record rather than recollection.

## Honest data scope
The analogue engine matches **only on periods present in the HDB** (crude ≈ 2018-12 → 2026). Each episode
below has a **computed dominant fingerprint** (modal `curveState`/`volState`/`liqState`/`rollPhase`/
`seasonPhase` + mean `slopePct`/`volPct`/`volumePct`) over its date range. The out-of-data lessons in the
last section are **narrative only — no data, never fingerprinted or matched.**

The episode date ranges + one-line risk-memory summaries are the spec in `.cfg.regime.episodes`; the full
failure-mode writeups live here, keyed by `driversKey`.

---

## In-data episodes (computed fingerprints; matchable)

### `crude_2020_covid` — COVID demand collapse + storage saturation
- **Range:** 2020-02-20 → 2020-06-30 · **driversKey:** `covid2020`
- **Drivers:** a global demand shock collapsed flat price while production was slow to cut; onshore/floating
  storage saturated, forcing a violent super-contango. WTI's May-2020 contract **printed negative
  (-$37) on 2020-04-20** as longs with no storage were forced out at expiry.
- **RISK MEMORY (what killed strategies here):**
  - **Long flat-price** carriers were run over by the crash and then by negative settlement mechanics.
  - **Short-vol / premium-selling** blew up as realized vol exploded.
  - **Calendar-spread / carry** books that were short the front dislocated as contango blew out far beyond
    any "fair" cost-of-carry — storage was the binding constraint, not arbitrage.
  - **Liquidity vanished** into the negative print; stops did not fill where expected. Position-sizing that
    assumed continuous markets was wrong.
  - Lesson: in a storage-saturation regime, the curve stops respecting cost-of-carry and gap risk dominates.

### `crude_2022_energyShock` — Russia/Ukraine supply shock
- **Range:** 2022-02-24 → 2022-08-31 · **driversKey:** `energyShock2022`
- **Drivers:** the invasion of Ukraine + sanctions risk drove a supply-side spike; the curve went sharply
  **backwardated** and realized + implied vol spiked.
- **RISK MEMORY:**
  - **Short-vol** and **mean-reversion shorts** (fading the rally) were run over by a one-way supply-driven
    move.
  - **Late trend/momentum chasers** who piled in near the peak **round-tripped** the reversal into H2-2022.
  - Backwardation roll-yield favoured longs — but the entry timing, not the carry, was the killer.
  - Lesson: a supply-shock backwardation rewards early trend and punishes both vol-sellers and late entrants.

### `crude_2023_rangebound` — OPEC+ cut-defended range
- **Range:** 2023-05-01 → 2023-12-31 · **driversKey:** `opecCuts2023`
- **Drivers:** OPEC+ supply cuts defended a price floor while soft demand capped the upside — a choppy,
  range-bound regime with no durable trend.
- **RISK MEMORY:**
  - **Trend / momentum** strategies **whipsawed** — breakouts repeatedly faded back into the range.
  - Only **modest carry** paid; directional conviction was punished in both directions.
  - Lesson: in a cut-defended range, fade breakouts and size trend-following down; the edge is mean-reversion
    inside the band, not the breakout.

---

## Out-of-data lessons (NARRATIVE ONLY — no data, not fingerprinted, not matched)

> These pre-date the HDB (crude starts ≈ 2018-12). They are **written memory only** — deliberately NOT
> entered as episodes with computed fingerprints, because there is no in-sample data to ground them. They
> are here so the lessons are not lost, not so the analogue engine can pretend to match them.

- **2008 — GFC demand collapse.** Flat price fell ~$147 → ~$30 in months; a deep-contango storage play
  emerged (the famous floating-storage trade). Killed: leveraged longs, and anyone short the contango too
  early before storage economics turned.
- **2014–2016 — shale supply glut.** A slow, grinding oversupply bear market into the high-$20s. Killed:
  premature bottom-pickers and carry books that underestimated how long contango could persist; rewarded
  patient short-carry and trend-following the downtrend.

*If/when earlier data is ingested, promote these to computed `regimeEpisodes` rather than leaving them as
narrative.*
