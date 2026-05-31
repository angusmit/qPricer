---
name: q-pricing
description: Use this skill when working on the user's q/kdb+ option pricing, commodity pricing, FDM, Greeks, scenario risk, validation, market data, AI-assisted workflows, maths testing, or quant technology project.
---

# q/kdb+ Option Pricing and Risk Framework Skill

You are helping with a q/kdb+ option pricing and risk analytics project.

The project goal is to build a front-office style option pricing and risk framework, not just a toy formula script.

The project starts with equity-style vanilla/FDM models as technical validation cases, but the longer-term target is a broader deterministic pricing and risk framework for commodity options, especially oil and power/electricity. Future extensions may include Merton jump diffusion, mean reversion models, mean-reverting jump models, stochastic volatility, and exotic options.

Core principles:
- q/kdb+ is the deterministic pricing, risk, data, and analytics engine.
- AI should assist with explanation, validation, debugging, code review, test generation, market data quality checks, maths knowledge testing, and workflow orchestration.
- Do not claim that AI produces the price. The pricing engine produces the price.
- Prefer small, testable changes over large rewrites.
- Preserve existing working code unless there is a clear bug or design problem.
- Always explain why a change is needed.
- Do not silently change pricing logic without adding or recommending tests.
- Treat AAPL/equity option data mainly as a technical test case, not the final target asset class.
- Keep the long-term roadmap aligned with commodity pricing and risk, especially oil and power/electricity.

## Project areas

The project may include:

- Black-Scholes pricing
- Binomial trees
- Explicit finite difference method
- Implicit finite difference method
- Crank-Nicolson finite difference method
- European options
- American options
- American option early exercise
- Greeks
- Scenario risk
- PnL explain
- Market data validation
- Option chain cleaning
- Volatility surface checks
- Calibration diagnostics
- q/kdb+ API layer
- AI explanation layer
- AI-assisted validation and debugging
- Commodity option pricing
- Oil option pricing
- Power/electricity option pricing
- Futures and forwards
- Convenience yield assumptions
- Seasonality
- Mean reversion
- Mean-reverting jump models
- Merton jump diffusion
- Jump risk diagnostics
- Stochastic volatility as a future extension
- Barrier options
- Asian options
- Digital options
- Basket options
- Lookback options
- Other exotic options as future extensions
- Quant maths learning and testing
- New technology testing and application
- Agentic workflows using tools such as MCP, Skills, coding agents, and AI assistants

## Long-term roadmap

The framework should develop in stages.

Stage 1: deterministic core pricing engine
- Black-Scholes baseline
- binomial tree
- FDM solvers
- Greeks
- scenario risk
- validation against known benchmarks

Stage 2: validation and market data layer
- option chain cleaning
- bid/ask validation
- stale quote detection
- implied volatility checks
- convergence tests
- benchmark comparisons
- risk report generation

Stage 3: commodity pricing extension
- oil options
- power/electricity options
- forward curve based pricing
- seasonality
- mean reversion
- jumps and spikes
- commodity-specific risk scenarios

Stage 4: advanced model extension
- Merton jump diffusion
- mean reversion
- mean-reverting jump diffusion
- stochastic volatility
- local volatility
- calibration diagnostics

Stage 5: exotic option extension
- Asian options
- barrier options
- digital options
- basket options
- lookback options
- path-dependent payoff framework

Stage 6: AI-assisted workflow layer
- natural language pricing request parser
- q error explainer
- validation assistant
- maths knowledge tester
- model selection assistant
- code review assistant
- market data quality explainer
- agent/MCP integration to call pricing tools
- README, documentation, and interview explanation generator

## New technology testing and application

When the user asks about new AI or software technology, evaluate it by how it can practically improve this q/kdb+ pricing project.

Relevant technologies may include:

- Claude Skills
- ChatGPT Projects / custom instructions
- GitHub Copilot
- Claude Code
- Codex or coding agents
- MCP
- agent workflows
- local tool calling
- Python wrappers around q
- q IPC
- Streamlit or dashboard UI
- automated tests
- CI/CD
- Docker
- local databases
- cloud deployment
- notebooks for research
- model validation reports

For each new technology, explain:

1. What it is.
2. Whether it is useful for this project.
3. Where it fits in the architecture.
4. What concrete feature it enables.
5. Whether it is needed now or later.
6. What the simplest first implementation would be.
7. What risks or distractions to avoid.

Example:
"MCP is useful later when you want an AI agent to call the q pricing engine directly. It is not needed before the q pricing functions and tests are stable."

## Maths knowledge testing and application

When helping with maths knowledge, do not only explain theory. Connect it to the pricing framework.

For each maths topic, include:

1. Core concept.
2. Why it matters in pricing or risk.
3. How it appears in the project.
4. What code module it relates to.
5. A simple test case.
6. A real desk-style interpretation.

Relevant maths and quant topics include:

- stochastic processes
- Brownian motion
- geometric Brownian motion
- Ito's lemma
- risk-neutral pricing
- martingales
- PDEs
- finite difference methods
- numerical stability
- interpolation
- root finding
- calibration
- optimisation
- Monte Carlo simulation
- variance reduction
- Greeks
- convexity
- volatility smile and skew
- local volatility
- jump diffusion
- mean reversion
- Ornstein-Uhlenbeck processes
- commodity forwards and futures
- convenience yield
- seasonality
- correlation
- copulas as a future topic
- credit-style hazard rate concepts if useful
- exotic option payoff logic

When the user asks to test their maths knowledge, use this format:

1. Quick concept check.
2. Small numerical example.
3. Link to pricing application.
4. q/Python implementation idea.
5. Interview-style explanation.

Example:
"If you understand why an American put should be worth at least the European put, you can test whether the FDM early exercise constraint is implemented correctly."

## q/kdb+ coding rules

When reviewing or writing q code:

1. Check atom vs list issues.
2. Check column type consistency.
3. Check table schema before applying functions.
4. Avoid unnecessary global state.
5. Avoid hidden mutation unless clearly documented.
6. Keep namespaces clear.
7. Prefer small functions with clear inputs and outputs.
8. Be careful with vector operations and implicit iteration.
9. Be careful with joins, keyed tables, and column ordering.
10. Do not over-engineer generic abstractions unless the repeated pattern is real.
11. Preserve working behaviour unless there is a clear reason to change it.
12. Suggest minimal patches first before proposing a full rewrite.
13. For every non-trivial code change, suggest at least one test.
14. Keep q code readable enough for portfolio review and interview discussion.

Common q errors to diagnose:
- `type`
- `length`
- `rank`
- `domain`
- `parse`
- `value`
- `nyi`

For each q error, explain:
- likely cause
- minimal reproducible example
- where to inspect
- suggested fix
- test to confirm

### q literal and comment syntax cautions

These cautions catch silently invalid q source that parses or runs without an obvious error and only manifests when a function is later called or proves to be missing. This list is **append-only and grows as new traps are found** — every milestone folds the lessons it hit back in here, so each trap is paid for once and never again.

**1. Typed numeric vector literal rule**

When creating typed numeric vectors in q, do not write invalid per-element suffix lists like `1 2 3f 4f 5f` or `1 2 3i 4i 5i`.

Use valid q syntax:

- `1 2 3 4 5f`
- `1 2 3 4 5i`
- or explicit casts like `"f"$1 2 3 4 5` and `"i"$1 2 3 4 5`.

Bad examples:

```q
1f 2f 3f
1i 2i 3i
1 2 3f 4f 5f
1 2 3i 4i 5i
0.25 0.5 1f 2f 5f
```

Good examples:

```q
1 2 3f
1 2 3i
1 2 3 4 5f
0.25 0.5 1 2 5f
"f"$1 2 3 4 5
```

Code review check: reject generated q code containing patterns like `1f 2f 3f`, `1i 2i 3i`, `1 2 3f 4f 5f`, or `1 2 3i 4i 5i`.

**2. q comment syntax caution**

A bare `/` line can start a multiline comment block until a standalone `\` is found, which may silently comment out following definitions. The file loads without error, but the affected functions never get defined and only fail when something tries to call or reference them.

Avoid bare `/` separator lines in q source files. Use `/ text` comments on every comment line instead.

Bad example:

```q
/
Helper functions
/
.my.ns.func:{[x] x+1}
```

In the bad example, the bare `/` on the first line opens a multiline comment block. Every subsequent line, including the `.my.ns.func` definition, is consumed by that block until a standalone `\` line appears.

Good example:

```q
/ Helper functions
.my.ns.func:{[x] x+1}
```

Code review check: scan q source files for any line consisting solely of `/` (no following text). If found, replace with `/ text` or remove the line.

**3. q is right-associative; parenthesise mixed-arity expressions**

q evaluates expressions right-to-left with uniform precedence. `2*metricCount=count rows` parses as `2*(metricCount=count rows)`, not `(2*metricCount)=count rows`. Likewise `string foo,")"` parses as `string (foo,")")` rather than `(string foo),")"` — wrap each conversion in parens, `(string foo),")"`, in every print-line concatenation (a bare `string x,"y"` threw `'type` and crashed the R12 demo's print loop). And a second difference (butterfly/curvature) must be written `(near+far)-2*mid`, NOT `near-2*mid+far` — the latter parses as `near-(2*mid+far)` and bit the R10 curve engine's curvature. Always bracket explicitly when mixing comparison, arithmetic, or function-application operators on the same line, especially in test assertions and print-line concatenations.

**4. `~` (match) is strict on type; `=` is not**

`1+til 4` produces a long vector `1 2 3 4`. `1 2 3 4i` is an int vector. `(1+til 4)~1 2 3 4i` is `0b` because match compares both value AND type. For test assertions on integer ranks/counts, prefer `(1+til n) ~ expectedLongVec` or use `=` if you only care about values. Same trap shows up for type codes: boolean atom is `-1h`, long atom is `-7h`, int atom is `-6h` — write the literal type code that matches your value's actual type.

**5. q lambdas do not capture outer locals like Python closures**

Inside a nested lambda, q can see its own args, its own locals, and globals — but NOT locals from an enclosing function body. If you define `requiredKeys` in an outer lambda and reference it inside a nested lambda, you'll get a `requiredKeys` error at call time. Pass the value in as an explicit argument and project it: `validateOne[;;requiredKeys]'[xs;ys]`.

**6. `select ... by ... from t` returns a keyed table; unkey before column access**

The result of `select … by key from t` is a keyed table. Direct column indexing (`kt`col`) on the key column can return unexpected shapes and surface as `type` errors in downstream `=` comparisons. Apply `0!` to convert to a plain table when callers expect ordinary column-name access — especially for test code or `exec … by` aggregations that flow into further `where`/`first` queries.

**7. Avoid q built-ins as variable/column/namespace names**

In addition to obvious cases (`value`, `count`, `select`, `key`, `type`, `string`, `first`, `last`, `sum`, `avg`, `min`, `max`, `exp`, `log`, `sqrt`, `asof`), the keyword `rank` is also reserved and will trigger `'assign` errors when used as a table-literal column name. Prefer descriptive names like `offenderRank`, `breachRank`, `optionValue`, `tradeCount`. Function parameters named `value` can also surface as `'match` errors at the lambda-definition site — rename to `atomValue` etc.

**8. q lambdas allow at most 8 explicit parameters**

`{[a;b;c;d;e;f;g;h;i] ...}` fails with `'params`. When you genuinely need more arguments — typically a wrapper that threads metadata plus existing state through several layers — bundle related args into dicts: pass `runMetadata:`runId`runDate`portfolioName!(…)` and `existingHistories:`limitHistory`summaryHistory!(…)` instead of nine bare positional params.

**9. `update col:scalar from t` on a single atom can misalign for strings**

For symbol/numeric/date atoms, q's `update col:atom from t` broadcasts the atom across rows correctly. For a string atom (char list), the same form treats it as a list and may error or assign char-by-char. The safe broadcast idiom is `count[t]#enlist atom` — this works uniformly for symbols, dates, ints, floats, and strings.

**10. `update col:atom from t` appends new columns at the end**

If the consumer cares about column order (tests, downstream `first cols tbl`, dashboard display), apply `xcols` to reorder after the `update`: `requiredColOrder xcols (update newCol:val from tbl)`.

**11. `n#x` cycles when n > count x**

`(-3)#enlist row` on a 1-row table returns a 3-row table by duplicating the single row, NOT an empty take and NOT an error. When you need a "take up to N rows, otherwise empty" semantics (e.g. computing a trend over the *prior* runs but the limit only has one observation), explicitly bound the take: `takeCount:lookbackRuns&priorAvailable; window:$[takeCount<=0;();(neg takeCount)#priorRows]`.

**12. Bare `/` lines: `\l` vs `value src` diverge**

A bare `/` line opens a multiline comment block that runs until the next bare `\` or EOF. The trap is that `value (read0 file)` does NOT honour multiline comments — it just evaluates statements — so a test that *evaluates* the file (via the project's `run_all_tests.q` harness) will pass, while *directly* running it with `q file.q` (which goes through `\l`) silently comments out everything after the bare `/` and the script returns exit 0 with no PASS line and no assertion failure. Always replace bare `/` with `/ text` or `/ -----`, especially in test headers.

**13. Seeded scan `f\[seed;list]` does NOT include the seed in the output**

`(+)\[10; 1 2 3]` returns `11 13 16` — three states from three list elements (cumulative sum from seed: 10+1=11, 11+2=13, 13+3=16). The seed `10` is the accumulator but is not part of the result. When you need to keep the initial state alongside the per-step results (typical for a backtest where step 0 IS the entry state), prepend explicitly: `states: enlist[initialState], f\[initialState; remainingRows]`. Or do `states: f\[initialState; pathRows]` and treat `pathRows[0]` as one of the step rows rather than the seed.

**14. `.[fn;arglist;handler]` needs a LIST of args, not a dict**

`.[fn;args;handler]` applies `fn` to `args` as a list (so `fn` is called with `args[0], args[1], ...`). If `fn` takes a single dict argument, pass `enlist dictArg`, not `dictArg` — otherwise q sees the dict's *values* as the arglist. Same for unary fn: use `enlist atomArg`.

**15. A list of dicts with identical keys auto-promotes to a TABLE**

`f\[seed;list]` returning a homogeneous list of state dicts is silently converted to a table by q. If you then prepend a single dict with a DIFFERENT key set via `enlist[d], that` you get `'mismatch` (table column-count vs the extra dict). Two safe fixes: (a) make the initial state share the post-step state shape exactly, so all states have the same keys — then prepending merges cleanly; (b) use a generic-list constructor like `(enlist d),`/`raze (enlist d;list)` only when you know q will keep the list-of-dicts representation. Option (a) is the more idiomatic and robust choice; it's also a useful design discipline for backtest state.

**16. `asc` sorts; the input order is not preserved**

`asc \`front\`back` returns `` `back`front `` (symbols sorted alphabetically), not the entry-order tuple. Test assertions that compare against the input-order tuple should compare against the SORTED expected tuple, or use `~` with the sorted form. For "every expected element is present" use `all expected in actual` instead. Same trap for numeric vectors: `asc 3 1 2` is `1 2 3`, not your input.

**17. Use a q TABLE as state for a mutable position book, not a list of dicts**

When a strategy's position evolves mid-path (legs created/retired), keep the legs as a TABLE column in state. Then qSQL handles the lifecycle naturally: `select … from legs where remainingTime <= rollThr` finds expiring legs; `select … from legs where not isExpiring` keeps survivors; concatenation `,` joins new legs; `update currentMark:unitPrice from legs` updates marks in bulk. Pricing each leg per step is a single `each` over the leg rows (giving a list of mark dicts), then `lj` on `legId xkey markedRows` brings the new marks alongside the leg specs. Avoid index-loops over legs; the table-based pattern composes with the rest of the qSQL the rest of the strategy already uses.

**18. `ss` is a q built-in (string-search); pick non-shadowing loop locals**

`ss` triggers a load-time `'ss` error when used as a local variable inside a `while` body. Same for other two-letter builtins (`ll`, `ml`, `mm`, **`sv`**, `vs`, etc.). For loop counters and strategy specs prefer descriptive names: `stratIdx`, `pathIdx`, `stratSpec`, `pathBundle`. SKILL note #7 covers the rule; this is a concrete trap caught while writing the ensemble runner. **`sv` specifically** (scalar-from-vector / join) bit a v0.52 test lambda whose params were `[d;x0;tenors;kappa;y0;sv;lv;cv]` — calling it raised `'nyi` (q resolved `sv` to the builtin). Capital-letter variants (`sV`, `lV`) are safe since names are case-sensitive; prefer descriptive `shortVolP`/`longVolP` regardless. Also note `csv` is a built-in CONSTANT (`csv` == `","`): assigning `csv:"..."` raises `'assign` (v0.54 trap) — use `csvText`.

**19. Initialise a "list of dicts" with `enlist d`, not `d`**

`synthRows:dict; synthRows,:enlist dict2` silently mis-amends the dict with the second dict's content because q's `,:` treats the lhs as a dict (key-merge), not a list. The correct seed for accumulating dicts is `synthRows:enlist dict0; synthRows,:enlist dict1` — both sides are lists of dicts and `,` concatenates them. This bit a test fixture where the assembled "rows table" came out as a single corrupted dict.

**20. `\`long$3.8` rounds half-up; use `floor` for truncating percentile indices**

`\`long$3.8` returns `4`, not `3` — q casts to long with rounding, not truncation. For a percentile-by-index formula like `idx: floor (pct%100) * (count vec - 1)`, use `floor` explicitly. Using `\`long$` made `p95` of `(1; 2; -1; 3; 0.5)` return `3` (the max) instead of `2` (the rank-3 value sorted). `floor` (or `_` for floor of a long) gives the conventional nearest-rank percentile.

**21. Inner lambdas can't see outer locals — inline or project tightly**

When you build helpers inside a function (`corOne:{...}; rowBuilder:{[i;pnls] corOne[pnls i;] each pnls}` inside a parent lambda), the inner `rowBuilder` will fail with `'corOne` at call time because q lambdas don't capture outer locals (rule #5). Either inline the helper into the place that needs it (`{[i;...] {...}[...] each ...}[;captured] each ...`), or pass it explicitly in via a projection. Most natural for matrix-builder-style code: write the whole double `each` inline with projection-captured arguments.

**22. `asof` is a built-in keyword, and a param sharing a column's name self-compares to always-true**

Two related parameter-naming traps in date/commodity code (both repeatedly hit across R9–R12): (a) `asof` is the q **as-of-join** keyword — never name a parameter or local `asof`; use **`asOf`** (the R9–R12 accessor/door/replay all standardised on `asOf`). (b) Never give a parameter the SAME name as a table column it filters on inside qSQL: a function `{[commodity] select from futures where commodity=commodity}` makes `where commodity=commodity` compare the COLUMN to itself — **always true**, silently returning every row. Use a distinct param name (`comm`, not `commodity`); the same applies to any `where col=col` where `col` is also the param. These parse and run without error — the bug is a wrong (over-broad) result set.

**23. A value immediately followed by `` `.ns.name `` parses the `.` as dot-apply, not a namespace path**

`(c\`price) .cfg.regime\`flatThreshold` does NOT index `c\`price` and then read a config value — q reads the ` .` as the **Apply/index operator** and evaluates `.[c\`price; cfg.regime\`flatThreshold]` (apply the price vector at that index), usually surfacing as `'type`/`'index` or a silent wrong value. Bind the namespaced config to a local FIRST (`thr:.cfg.regime\`flatThreshold; ... use thr`), or fully parenthesise both operands. This bit the R10 curve engine reading `.cfg.regime` thresholds.

**24. qSQL (`select`/`exec`/`update`) on a LOCAL table inside a lambda throws `'assign`**

`{[run] steps:run\`steps; fills:select from steps where filledQty<>0f; …}` raises `'assign` at PARSE time: q's qSQL templates resolve the table name against the GLOBAL namespace, not a function-local, so a local table is not a valid target. Two fixes: (a) plain vector ops — `mask:0f<>steps\`filledQty; fillQ:(steps\`activeContract) where mask` (filter columns, not the table); or (b) the functional form `?[steps;enlist(<>;\`filledQty;0f);0b;()]`. This bit the R13 evidence audit (and its test harness) repeatedly — prefer column/`where`-mask extraction over `select … from <local>` inside any helper lambda. (Also recurring: a backtick-INDEX followed by a join, `dict\`reason,")"`, parses as `dict[\`reason,")"]` — the `,` binds the symbol to the string BEFORE the index; parenthesise the index: `(dict\`reason),")"`.)

**25. The `.mm` / `.dd` / `.month` dot-accessors are ATOM-ONLY — they error on a temporal VECTOR**

`someDate.mm` gives the month 1-12 for a date ATOM, but `dateVec.mm` raises `'.mm` on a date VECTOR. For a vector calendar-month (e.g. grouping by same calendar month), cast instead: `(\`month$dateVec) mod 12` (0-11, Jan=0 — equality is all that matters for same-month grouping). This bit the R14 seasonality layer building the same-calendar-month series. (Same atom-only caveat applies to `.dd`/`.month`/`.year` on vectors — prefer the `\`month$`/`\`date$` casts or `mod` arithmetic for vectorised temporal-part extraction.)

**26. The test runner `value`s every test in ONE process — globals PERSIST across tests; isolate state a test depends on**

`tests/run_all_tests.q` evaluates each test file in the SAME q process, so any GLOBAL a test sets leaks into every later test (in suite order). A test that relies on a LIBRARY default of a mutable global must RESET it at the top, or a prior group's leftovers silently flip its result (passes standalone, fails in the full suite — the worst kind). The repeat offenders: `.gov.hypoTbl`/`.gov.trialTbl` (reset via `.gov.hypoTbl:.gov.__emptyHypotheses[]` — and `.gov.register`'s upsert can `'type` on a type-incompatible leftover row, R13/R16), the cards store `modelCards` (the cards tests leave a synthetic one; restore the real store with `modelCards:.cfg.cards` so `.cards.gateReady` finds the real cards, R16), and `futures` (every commodity/replay/season test sets its own — assume nothing about it on entry). Also use UNIQUE hypoIds when registering. Rule: if a test reads a mutable global it did not itself set this run, reset it first.

## q/kdb+ naming and namespace rules

When reviewing or writing q code, enforce clear naming and namespace hygiene.

### Avoid reserved words and built-in q names

Do not create variables, functions, namespaces, or table columns that shadow q keywords, q operators, system functions, or common built-in names.

Avoid names such as:

- `select`
- `exec`
- `update`
- `delete`
- `from`
- `where`
- `by`
- `x`
- `y`
- `z`
- `i`
- `j`
- `k`
- `key`
- `value`
- `count`
- `first`
- `last`
- `sum`
- `avg`
- `max`
- `min`
- `prd`
- `var`
- `dev`
- `sqrt`
- `exp`
- `log`
- `til`
- `raze`
- `enlist`
- `flip`
- `cols`
- `meta`
- `type`
- `rank`
- `string`
- `symbol`
- `date`
- `time`
- `each`
- `peach`
- `prior`
- `next`
- `prev`
- `over`
- `scan`
- `within`
- `like`
- `in`
- `not`
- `null`
- `tables`
- `views`
- `system`
- `hopen`
- `hsym`
- `read0`
- `read1`
- `save`
- `load`
- `.Q`
- `.q`
- `.z`
- `.h`
- `.j`
- `.o`

If a name conflicts with a q keyword, built-in function, or common system namespace, rename it to a more precise project-specific name.

Bad examples:

```q
value:100
count:{x+y}
select:{x}
price:{x}
d:{x+y}
```

Better examples:

```q
optionValue:100f
tradeCount:count tradeTable
calculateOptionPrice:{[pricingConfig] ... }
discountFactor:{[riskFreeRate;timeToExpiry] exp neg riskFreeRate*timeToExpiry}
```

### Namespace rules

Use explicit project namespaces.

Good namespace examples:

```q
.pricing.blackScholes.priceEuropeanCall
.pricing.fdm.priceAmericanPut
.pricing.greeks.calculateDelta
.pricing.validation.validateOptionConfig
.pricing.market.cleanOptionChain
.pricing.risk.calculateScenarioRisk
.pricing.commodity.priceOilOption
.pricing.commodity.pricePowerOption
.pricing.jump.priceMertonOption
.pricing.meanReversion.priceMeanRevertingOption
.pricing.exotic.priceAsianOption
```

Avoid unclear or overly short namespaces:

```q
.d.price
.p.f
.x.calc
.u.test
```

Do not use `.d` as a project namespace because it is too vague. Prefer descriptive namespaces such as:

```q
.pricing
.market
.validation
.greeks
.risk
.solver
.model
.product
.commodity
.jump
.meanReversion
.exotic
.test
```

### Public and private function convention

Public functions should be placed near the top of the file.

Private helper functions should use double underscore after the namespace:

```q
.pricing.fdm.priceAmericanPut:{[pricingConfig] ... }
.pricing.fdm.__buildGrid:{[gridConfig] ... }
.pricing.fdm.__applyBoundaryConditions:{[grid;pricingConfig] ... }
.pricing.fdm.__stepBackward:{[grid;modelConfig] ... }
```

Public functions should represent the external API of the module.

Private functions should only support public functions in the same namespace.

When reviewing code, check whether private helper functions are accidentally being used as public APIs.

### Naming style

Use descriptive camelCase names for variables, functions, config keys, and table columns.

Good examples:

```q
spotPrice
strikePrice
timeToExpiry
riskFreeRate
dividendYield
volatility
optionType
exerciseStyle
gridSize
timeStepCount
priceEuropeanOption
calculateScenarioRisk
validatePricingConfig
meanReversionSpeed
longRunMean
jumpIntensity
jumpMean
jumpVolatility
commodityForwardPrice
seasonalityFactor
```

Avoid short or unclear names except in very small local mathematical expressions.

Bad examples:

```q
s
k
t
r
q
v
n
nj
cp
res
tmp
foo
bar
```

Better examples:

```q
spotPrice
strikePrice
timeToExpiry
riskFreeRate
dividendYield
volatility
spaceStepCount
timeStepCount
callPutFlag
pricingResult
temporaryGrid
```

Short names such as `i`, `j`, `k`, `x`, and `y` are acceptable only for tiny local loops or anonymous functions, but not for public API, table columns, config keys, or important model variables.

### Table column naming

Use precise camelCase column names.

Good examples:

```q
tradeId
symbol
optionType
exerciseStyle
spotPrice
strikePrice
timeToExpiry
riskFreeRate
dividendYield
volatility
modelPrice
marketMidPrice
delta
gamma
theta
vega
rho
scenarioName
scenarioPnL
validationStatus
validationReason
commodityName
deliveryMonth
forwardPrice
seasonalityFactor
jumpIntensity
meanReversionSpeed
```

Avoid vague columns:

```q
s
k
t
v
px
res
flag
x
y
```

### Code review behaviour for naming

When reviewing code, always check:

1. whether any custom name shadows a q built-in;
2. whether namespaces are explicit and meaningful;
3. whether public functions are easy to find;
4. whether private helper functions use `.__`;
5. whether variable names are precise enough for pricing and risk code;
6. whether short names make the code harder to review;
7. whether table columns are clear enough for downstream analytics;
8. whether config keys are understandable without reading the implementation.

If names are unclear, suggest a rename table like:

| Current name | Suggested name | Reason |
|---|---|---|
| `s` | `spotPrice` | More explicit |
| `k` | `strikePrice` | Avoid unclear mathematical shorthand |
| `v` | `volatility` or `optionValue` | Ambiguous |
| `.d.price` | `.pricing.fdm.priceOption` | Explicit namespace and purpose |

## FDM validation checklist

When working on FDM pricing:

1. Confirm product type:
    - European call
    - European put
    - American call
    - American put

2. Confirm model assumptions:
    - spot
    - strike
    - maturity
    - volatility
    - rate
    - dividend yield
    - boundary conditions
    - grid size
    - time steps

3. Check numerical method:
    - explicit
    - implicit
    - Crank-Nicolson

4. Check boundary conditions:
    - call lower boundary near zero
    - call upper boundary approaches S - discounted strike
    - put lower boundary approaches discounted strike
    - put upper boundary near zero

5. For European options:
    - compare with Black-Scholes
    - check convergence as grid becomes finer
    - check put-call parity where applicable

6. For American options:
    - American put should be worth at least European put
    - early exercise premium should be non-negative
    - apply max(option value, intrinsic value) at each time step
    - compare with binomial benchmark if available

7. Check stability:
    - explicit method has stability constraints
    - implicit and Crank-Nicolson should be more stable
    - bad grid choices can cause oscillation or negative prices

8. Check output reasonableness:
    - option price should not be negative
    - call price should not exceed spot price
    - put price should not exceed discounted strike too aggressively unless American exercise and assumptions justify it
    - price should move in the expected direction when spot, volatility, or maturity changes

## Commodity pricing checklist

When working on commodity pricing, check:

1. Whether the underlying is spot, futures, or forward based.
2. Whether the option is written on spot, futures, average price, spread, or delivery period.
3. Whether the asset is oil, gas, power/electricity, metals, or another commodity.
4. Whether the model should include:
    - forward curve
    - convenience yield
    - storage costs
    - seasonality
    - mean reversion
    - jumps or spikes
    - delivery period effects
5. For oil:
    - consider futures curve and convenience yield assumptions.
    - consider jump risk for supply shocks.
6. For power/electricity:
    - remember power is non-storable.
    - consider spikes, seasonality, and mean reversion.
    - GBM may be only a baseline, not a realistic final model.

When explaining commodity model choice, link the model to the asset behaviour.

Example:
"Mean reversion is more relevant for power than pure GBM because electricity prices can spike and then revert as supply/demand normalises."

## Jump diffusion and mean reversion checklist

When working on Merton jump diffusion:

- identify jump intensity
- identify jump size mean
- identify jump size volatility
- check how jumps affect return distribution
- compare against Black-Scholes baseline
- test limiting case where jump intensity is zero
- explain whether jumps help capture smile/skew or tail risk

When working on mean reversion:

- identify long-run mean
- identify mean reversion speed
- identify volatility
- check whether the process is on price, log-price, spread, or factor
- test limiting cases
- explain why mean reversion is useful for commodities

When working on mean-reverting jump models:

- combine reversion and jumps carefully
- explain which part captures normalisation and which part captures spikes
- test that jumps do not break validation assumptions
- include clear parameter names and calibration diagnostics

## Exotic option checklist

When working on exotic options:

1. Identify payoff type:
    - Asian
    - barrier
    - digital
    - basket
    - lookback
    - spread
    - other path-dependent payoff

2. Identify path dependency:
    - terminal only
    - running average
    - barrier monitoring
    - max/min dependence
    - multiple underlyings

3. Identify suitable numerical method:
    - closed form if available
    - binomial/trinomial tree
    - FDM
    - Monte Carlo
    - PDE with extra state variables
    - approximation

4. Validation:
    - compare with vanilla special case where possible
    - test limiting cases
    - test monotonicity
    - test payoff shape
    - test path simulation logic

5. Explanation:
    - describe payoff clearly
    - explain why vanilla pricing is insufficient
    - explain key risk drivers
    - explain model limitations

## Pricing result explanation

When explaining a pricing result, include:

- model used
- option price
- whether result is reasonable
- comparison with benchmark if available
- main assumptions
- Greeks interpretation
- scenario risk interpretation
- limitations

Use front-office / desk-style language, but keep it understandable.

Example:
"The American put is above the European put, which is expected because early exercise has value. The positive early exercise premium suggests the solver is applying the exercise constraint correctly."

## Greeks explanation

When explaining Greeks:

- Delta: directional exposure to spot.
- Gamma: convexity of delta; explains how delta changes as spot moves.
- Theta: time decay; explains how option value changes as time passes.
- Vega: volatility exposure; explains sensitivity to implied volatility or model volatility.
- Rho: interest-rate sensitivity.

When explaining a portfolio or option result, translate Greeks into risk language.

Example:
"The option has negative delta, so it benefits when the underlying falls. Positive gamma means the position becomes more directionally short as spot falls and less short as spot rises. Positive vega means the option gains value when volatility increases."

For commodity options, also consider:
- sensitivity to forward curve shifts
- sensitivity to seasonality
- sensitivity to mean reversion parameters
- sensitivity to jump intensity
- sensitivity to delivery month

## Scenario risk explanation

When explaining scenario risk:

- Explain which scenario causes the largest gain.
- Explain which scenario causes the largest loss.
- Link scenario PnL to Greeks.
- Mention nonlinear effects where gamma, jumps, mean reversion, or early exercise may matter.
- Do not overstate precision. Scenario results depend on the model and assumptions.

Example:
"The spot-down scenario is positive because the position is a put with negative delta. The vol-up scenario is also positive because the option has positive vega. The spot-up scenario is negative because the put loses intrinsic and time value as the underlying rises."

For commodity options, include scenarios such as:
- parallel forward curve shift
- front-month shock
- back-month shock
- seasonality shock
- volatility shock
- jump intensity shock
- mean reversion speed shock

## Market data validation checklist

When reviewing option chain or market data, check:

- missing spot price
- missing strike
- missing expiry
- missing option type
- missing bid or ask
- bid greater than ask
- zero or negative ask
- negative bid
- stale quote
- very wide bid-ask spread
- missing implied volatility
- unrealistic implied volatility
- missing open interest or volume
- duplicate contracts
- inconsistent expiry format
- invalid call or put flag
- inconsistent moneyness
- suspicious volatility surface points
- calendar arbitrage issues
- vertical spread arbitrage issues where relevant

For commodity market data, also check:
- futures curve completeness
- missing delivery months
- inconsistent contract codes
- roll date handling
- seasonality pattern consistency
- negative or extreme prices where applicable
- spikes or outliers
- liquidity by contract month
- settlement vs intraday price differences

When explaining market data problems, separate:

- hard errors that must be removed;
- soft warnings that should be reviewed;
- modelling limitations that may require better assumptions.

## Volatility surface and calibration diagnostics

When reviewing volatility surface or calibration output:

1. Check whether the surface is smooth enough for the model.
2. Check whether there are obvious outliers.
3. Check whether calibration error is concentrated in specific strikes or expiries.
4. Explain whether error is likely due to:
    - poor liquidity
    - stale quotes
    - wide bid-ask spread
    - missing dividends
    - missing forward curve assumptions
    - flat-vol model limitation
    - local-vol sensitivity
    - jump risk
    - commodity seasonality
    - mean reversion
    - bad interpolation
5. Recommend whether to clean, exclude, smooth, or investigate points.

Do not treat every model-market difference as true mispricing.

## AI-assisted project features

Good AI features to suggest:

- natural language pricing request parser
- q error explainer
- FDM validation assistant
- commodity model selection assistant
- market data quality explainer
- Greeks and scenario risk explanation
- model selection assistant
- calibration diagnostics
- README and portfolio write-up generator
- test case generation assistant
- code review assistant for q/kdb+ naming, namespaces, and structure
- maths knowledge tester
- interview-style question generator
- MCP tool layer for calling q pricing functions
- agent workflow for price, validate, explain, and report

Avoid suggesting a black-box AI model as the primary pricing engine at the beginning.

The preferred project positioning is:
"q/kdb+ deterministic pricing and risk engine with an AI-assisted explanation, validation, maths testing, and workflow layer."

## Agent and MCP roadmap

If the user asks about agents or MCP, explain them in relation to this project.

MCP should be used later to expose stable q/Python functions as tools, such as:

- priceOption
- calculateGreeks
- runScenarioRisk
- validateOptionChain
- runFdmConvergenceTest
- calibrateModel
- explainQError
- generateRiskReport

An AI agent can then orchestrate:

1. parse user request
2. validate inputs
3. call pricing engine
4. calculate Greeks
5. run scenarios
6. compare against benchmarks
7. explain result
8. suggest tests or next improvements

Do not introduce MCP before the core pricing APIs are stable enough to call.

## Code review response format

When reviewing code, structure the response like this:

1. Summary of what the code is trying to do.
2. Main issue or risk.
3. Specific q/kdb+ issue if any.
4. Suggested minimal patch.
5. Why the patch works.
6. Test to confirm.
7. Optional improvement if the user wants to refactor later.

If there are naming issues, include a rename table.

If there are pricing issues, include a validation checklist.

If there are q errors, include a minimal reproducible example where possible.

## Maths testing response format

When testing maths knowledge, structure the response like this:

1. Concept question.
2. Expected intuition.
3. Small calculation or derivation.
4. Link to pricing/risk.
5. q or Python implementation idea.
6. Interview-style answer.

Example:
"Explain why increasing volatility increases the value of a vanilla option. Then connect this to positive vega and a vol-up scenario."

## New technology evaluation response format

When evaluating new technology, structure the response like this:

1. What it is.
2. What problem it solves.
3. Whether it helps this project.
4. Where it fits:
    - coding
    - testing
    - pricing
    - validation
    - explanation
    - deployment
    - agent workflow
5. Use now or later.
6. Smallest useful implementation.
7. Risk or distraction.

## Do not do

Do not:
- rewrite the whole project without being asked;
- replace q pricing logic with Python unless specifically requested;
- hide changes behind vague explanations;
- invent validation results that were not run;
- claim the model is correct without benchmark or test evidence;
- create unclear namespaces like `.d`, `.p`, or `.x`;
- create public functions with one-letter names;
- use q built-in names as custom variables or functions;
- overcomplicate the design before the basic pricing and validation work;
- treat AAPL/equity options as the final target if the discussion is about the long-term project;
- recommend complex AI/agent infrastructure before the deterministic pricing and validation core is stable.
