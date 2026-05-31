/ templates/template.q - the problem-TEMPLATE abstraction (.template.*) (v0.70, Research OS R6)
/ ----------------------------------------------------------------------------
/ ARCHITECTURE.md Part II (the three plug-in kinds) / 13-R6. A problem template is a
/ research SHAPE that COMPOSES capabilities, declares its OWN inputs/outputs/validation,
/ and produces a per-day PnL series that flows into the universal gov gates (R3/R3b).
/ It is the THIRD plug-in kind (CAPABILITY = R2 contracts, KNOWLEDGE = R5 cards,
/ PROBLEM-TEMPLATE = here), registered through R2's `.registry.*` as a `template` kind, so
/ a template is a first-class, conformance-checked, card-able plug-in.
/ -
/ LAYER: templates/ composes signals//backtest//execution/ (loads AFTER them), is
/ INDEPENDENT of gov/ (no mutual import - a template PRODUCES a PnL series + a validation
/ manifest; the gov gates CONSUME a PnL series). Loads after backtest//portfolio/, before
/ gov/. Never opens the HDB at import.
/ -
/ A template's run returns `pnl`validation`meta:
/   pnl        - a (date;pnl) table: the per-day NET PnL series the gov gates will judge,
/   validation - the TEMPLATE-SPECIFIC check results (a directional template has none beyond
/                the universal gates; relativeValue adds a spread-stationarity gate),
/   meta       - shape + provenance.
/ The `template` contract enforces this uniform OUTPUT; inputs are template-specific
/ (declared in each manifest for documentation, not contract-required).
/ ----------------------------------------------------------------------------

/ The template-kind contract: the uniform OUTPUT shape (inputs vary by template).
.template.contract:`version`requiredIn`requiredOut!(
    1;
    ()!();
    `pnl`validation`meta!`table`dict`dict);
.registry.new[`template; .template.contract];

/ Per-kind registry (thin specialisation over .registry.*, the `template` kind).
.template.register:{[name;fn;manifest] .registry.register[`template;name;fn;manifest]};
.template.get:{[name] .registry.get[`template;name]};
.template.list:{[] .registry.list `template};
.template.conforms:{[name] .registry.conforms[`template;name]};

/ Run a registered template on its inputs -> `pnl`validation`meta.
.template.run:{[name;inputs] ((.registry.get[`template;name])`fn) inputs};

/ ── template #1: directionalSignal (the FAITHFUL WRAPPER - byte-identity proof) ──

/ Expresses the existing directional shape as a template by DELEGATING to the unchanged
/ signal + backtest functions - it does NOT reimplement the logic, so its per-day PnL is
/ identical BY CONSTRUCTION to the existing backtest path. inputs:
/   strategy (sym), trade (dict), stratCfg (dict), model (dict), optional fdmCfg (dict),
/   and EITHER a prebuilt `path (the strategy path table) OR `curveHistory + `sigCfg
/   (the commodity path is then built via the unchanged .strategy.path.commoditySignals).
/   Optional dateFrom/dateTo restrict the curve history (the gov runner seam).
.template.directional.run:{[inputs]
    path:$[`path in key inputs;
        inputs`path;
        [ch:inputs`curveHistory;
         if[(`dateFrom in key inputs) and `dateTo in key inputs;
             ch:select from ch where asofDate within (inputs`dateFrom;inputs`dateTo)];
         (.strategy.path.commoditySignals[ch; inputs`sigCfg])`path]];
    fdmCfg:$[`fdmCfg in key inputs; inputs`fdmCfg; ()!()];
    res:(.strategy.runAndSummarize[inputs`strategy; inputs`trade; path; inputs`model; fdmCfg; inputs`stratCfg])`result;
    ok:select from res where status=`OK;
    pnl:select date:stepDate, pnl:stepPnl from ok;
    `pnl`validation`meta!(
        pnl;
        (enlist `shapeGate)!enlist `passed`detail!(1b;"directionalSignal: no template-specific gate; the universal gov gates judge it");
        `shape`strategy`nDays`edgeSource!(`directionalSignal; inputs`strategy; count pnl; `riskPremium))
 };

.template.register[`directionalSignal; .template.directional.run;
    `contractVersion`description`shape`composes`in`out!(
        1;
        "directional signal strategy - a faithful wrapper delegating to the unchanged signal + backtest path";
        `directionalSignal;
        `commoditySignalPath;
        `strategy`trade`stratCfg!`symbol`dict`dict;
        `pnl`validation`meta!`table`dict`dict)];

-1 "template.q loaded - .template.* abstraction + directionalSignal (faithful wrapper) ready";
