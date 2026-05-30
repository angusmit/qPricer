---
name: readme-audit
description: Audit and refresh every README.md in the repo so each accurately and completely describes its directory's current contents (modules, namespaces, key API, dependencies). Use when READMEs have drifted or after a milestone that added or changed modules. Docs-only; touches no code.
---

# README audit — bring all READMEs to current state

DOCS ONLY. The only files changed must be *.md. Follow the standing milestone conventions in
CLAUDE.md (shell hygiene, scratch, no cleanup, no mkdir, no orphan q).

## Procedure
1. Enumerate every README: `git ls-files '*README.md'` (+ check for untracked), plus the root README.md.
2. For each, inspect the ACTUAL current contents of its directory — the .q modules present, the
   namespaces they define, the key public functions.
3. Rewrite each layer README to MATCH reality, uniform structure:
    - Purpose (1-2 sentences)
    - Dependencies (which layers it depends on; deps flow downward only, ARCHITECTURE.md §1)
    - Modules (every .q file, one line each)
    - Key API (the main namespaces/functions callers use)
    - Notes (gotchas / config keys / reserved-layer status)
      Root README.md: project overview, current layer map, Quick Start (`\l core/init.q`), accurate tree.
4. Cross-check each README against CLAUDE.md and ARCHITECTURE.md; where a README contradicts the
   current state, fix the README (do NOT silently edit CLAUDE.md/ARCHITECTURE.md — flag contradictions).

## Boundaries
- tests/ and all .q/.cfg files are OFF LIMITS — change *.md only. Don't reduce accuracy; complete and correct.

## Verify
- `git diff --name-only` shows ONLY *.md.
- Each README lists every .q module currently in its dir + the real API.
- A fresh `\l core/init.q </dev/null` still loads. No full-suite run needed (docs-only).

## Report
- Every README found and what was stale/missing vs updated (or "already current"); confirm the diff
  is *.md-only; any README-vs-top-doc contradictions; suggested docs-only commit message.