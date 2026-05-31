# Agent orchestration contract (Research OS R7)

Six **bounded** research-agent roles enact one disciplined loop on top of the qFDM Research OS (R1–R6 + the v0.71 wiring). The defining principle is **boundedness**: the agents propose, gate, challenge, record, curate, and coordinate — **none is an autonomous trader, none allocates capital, and the human holds every go/no-go**. qFDM is a **batch research / simulation system**: there is no live trading and no order execution; the agents produce research output (an honest verdict + the case) for a human.

## The roles
| Role | File | Does |
|---|---|---|
| Researcher | `researcher.md` | Proposes + pre-registers the hypothesis (thesis + edge source, a-priori); instantiates a template |
| Curator | `curator.md` | Ensures the capability is carded with a populated failure-mode field; maintains cards (R5) + the regime library (R4) |
| Validator | `validator.md` | Runs `.cards.gatedRun` → `.gov.runFull`; reports the honest verdict; the only role that touches the holdout |
| Skeptic | `skeptic.md` | Surfaces the regime risk memory; hunts overfitting / post-hoc / priced-in; finds why it won't work |
| Logger | `logger.md` | Guards the honest record — every trial logged (honest N), holdout spent ≤ once, no cherry-picking |
| Lead | `lead.md` | Orchestrates the loop, enforces the merge gate, packages the case, **escalates every go/no-go to the human** |

(Maps to ARCHITECTURE.md §11.7's conceptual stages: Rationale/Proposer → History/Regime → Model → Validation → Skeptic → Ledger, with the Lead as the coordinator/merge-gate enforcer.)

## The loop
```
researcher proposes (register, a-priori)
  -> curator ensures carded (failure modes documented; else REFUSED before any gate)
  -> validator gates           (.cards.gatedRun -> .gov.runFull: cost -> deflation -> walk-forward -> one-shot holdout)
  -> skeptic challenges         (riskMemory annotation: "resembles <episode>; here's what killed strategies then")
  -> logger records             (every trial -> honest N; holdout spent <= once)
  -> lead packages + ESCALATES to the human   (verdict + tradeable fact + DSR/holdout + risk memory + the case)
```
The single-process reference is `.workflow.run[proposal]` (`workflow/workflow.q`); the demo is `apps/examples/research_loop.q`.

## The hard invariants (each non-negotiable)
1. **Agents are bounded** — no autonomous trading, no autonomous capital allocation; the **human** holds every go/no-go.
2. **The sealed holdout** is touched only by the validator's one-shot gate (R3b) — no role peeks, ever; a second look returns the recorded verdict.
3. **The ledger logs EVERY trial** (the honest N); survivors are never cherry-picked.
4. **"Tradeable" requires passing ALL gates** (the R3b fail-safe) **AND** a human sign-off — passing the gates is necessary, never sufficient.
5. **A capability is gateable only when carded** with a populated failure-mode field (the v0.71 carded gating; `.cards.gatedRun` refuses otherwise).
6. **Capital staging is paper → small → scaled**, advancing only on **explicit human sign-off** at each stage. (This refers to the human's eventual real-world deployment decision, which the system informs but never performs.)
7. **The merge gate** for any agent-produced code change is the **full suite green + byte-identical canonicals**, run via the test-runner / main thread.

The discipline lives in the **spine**, not in the agents' goodwill: the holdout is technically sealed, the ledger is append-only, gating refuses undocumented capabilities, and nothing is marked tradeable without every gate — so scaling the workforce never loosens the discipline.

## Deferred (the next step)
The **multi-agent orchestration** (Claude Code worktrees / Agent Teams / Dynamic Workflows) is a named deferral. R7 ships the role definitions, this contract, and the single-process `.workflow.run` loop as the reference the later orchestration parallelises — **behind the same merge gate**, with `CONTRACTS.md` + `CLAUDE.md` + the q-pricing skill as the shared constitution.
