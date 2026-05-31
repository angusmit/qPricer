# workflow/ — The single-process research loop (ARCHITECTURE.md §11.7, R7)

## Purpose
A thin HIGH layer (`.workflow.*`) that COMPOSES the already-built functions into the bounded research loop the six agent roles (`agents/*.md`) enact: researcher proposes → curator ensures carded → validator gates → skeptic annotates → logger records → lead packages + escalates to the human. It is the single-process reference the later multi-agent orchestration will parallelise.

## Dependencies
COMPOSES `.gov.*` (R3/R3b) / `.cards.*` (R5 + the v0.71 carded gating) / `.template.*` (R6) / `.regime.*` (R4) — **downward only** (it is the top functional layer; loads after `cards/`, before `apps/`). It REUSES; it reimplements nothing, adds no compute path, and opens no HDB at import.

## Modules
- `workflow.q` — `.workflow.run` (the loop) + `.workflow.__packet` (the human-escalation packet).

## Key API
- `.workflow.run[proposal]` → a **human-escalation packet**. `proposal` is a dict: `hypoId`, `thesis`, `edgeSource`, `instruments`, `claimedRegimes`, `capName` (the capability to gate — must be carded + registered), `runner` (a `runner[from;to]→(date;pnl)` callback), optional `axis` (default `curveState`).
  - Steps: (1) researcher — `.gov.register` (pre-register, a-priori); (2) curator — `.cards.gateReady` (REFUSE here if the capability has no populated failure-mode field — gates NOT run, holdout untouched); (3–5) validator + skeptic + logger — `.cards.gatedRun` → `.gov.runFull` (the one-shot-holdout cascade, the `riskMemory` annotation, and the trial logging are already composed there); (6) lead — assemble the packet.
- The packet: `hypoId`, `capability`, `carded`, `gatesRun`, `anyTradeable` (a DERIVED fact — true iff some bucket passed ALL gates), `verdicts` (the per-bucket table), `riskMemory` (the skeptic annotation), `decision` (always `escalateToHuman`), `humanSignOffRequired` (always `1b`), `stage` (`paper`), `summary`.

## Bounded by construction
The loop has **no path** that allocates capital, trades, or marks something tradeable without all gates passing. `anyTradeable` is a fact, not a decision; the packet carries **no** `deploy`/`allocate`/`autoTrade` field and **always** defers to the human (`tests/workflow/test_workflow_bounded.q` proves both end-to-end). qFDM is a batch research system — there is no order execution.

## Notes
- ADDITIVE: composes existing entry points, changes no compute path, so the full suite stays byte-identical. Demo: `apps/examples/research_loop.q`.
- The multi-agent orchestration (worktrees / Agent Teams / Dynamic Workflows) is a named deferral — see `agents/ORCHESTRATION.md`; it parallelises this same loop behind the same merge gate.
