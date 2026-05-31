/ execution/execution.q - daily fill-and-cost execution simulation (.exec.*) (v0.60)
/ ----------------------------------------------------------------------------
/ Migration step 4 (ARCHITECTURE.md sections 5 / 9). Turns a strategy's TARGET
/ position change (an order) into a simulated FILL (with optional slippage and a
/ liquidity participation cap) and explicit COST components, so a backtest can
/ book REALIZED (net-of-execution) PnL and report GROSS vs NET performance.
/ -
/ Matched to the data: daily SETTLE prices only (no bid-ask, no intraday), so
/ this is a daily fill-and-cost layer, NOT a tick / limit-order-book simulator.
/ -
/ Units: positions here are the backtest's VOL-TARGETED notional scales (not
/ integer contracts), so the proportional cost is rate * |filled turnover| (no
/ price multiply) - exactly the legacy txnCostRate * turnover. Slippage IS a
/ price-scaled adverse cost (bps of traded notional = |filled| * refPrice * bps).
/ -
/ BYTE-IDENTICAL DEFAULT: with the default config (proportional cost at the
/ strategy's own txnCostRate, zero slippage, zero fixed, no participation cap)
/ the order fills fully and totalCost == legacy cost, so a wired backtest
/ reproduces the pre-execution-layer numbers to the digit. Realism is opt-in.
/ -
/ Generic + stateless: .exec.fill is a pure function the equity engine can also
/ adopt later (step 4b); this step wires only the commodity backtest.
/ ----------------------------------------------------------------------------

/ Default execution config (returns the .cfg.exec dict populated by config/base.q).
.exec.defaultConfig:{[] .cfg.exec};

/ Resolve the execution config for a strategy run: start from .cfg.exec, apply an
/ optional per-run `exec` sub-config from stratCfg, then force proportionalRate to
/ the strategy's own txnCostRate (the legacy proportional rate is authoritative,
/ so the default path reproduces the legacy cost regardless of .cfg.exec's default
/ rate). Realism levers (slippageBps / fixedPerTrade / participationCap / impact)
/ come from the per-run `exec` sub-config and are off by default.
.exec.__resolve:{[stratCfg]
    cfg:.cfg.exec;
    if[`exec in key stratCfg; cfg:cfg,stratCfg`exec];
    if[`txnCostRate in key stratCfg; cfg[`proportionalRate]:stratCfg`txnCostRate];
    cfg
 };

/ Simulate one fill. order is the SIGNED desired position change (target-current).
/ ctx carries: refPrice (the fill reference, = the price the return is booked on),
/ optional barVolume / adv (for the participation cap + impact) and volatility
/ (for vol-scaled slippage), and currentPos. cfg is a resolved .cfg.exec dict.
/ Returns a dict: filledQty (signed), fillPrice, proportionalCost, slippageCost,
/ fixedCost, totalCost. Pure / no state. Known-answer testable.
.exec.fill:{[order;ctx;cfg]
    refPrice:ctx`refPrice;
    side:signum order;
    absOrder:abs order;
    / participation cap: fill fully unless a finite cap AND a finite bar volume bind.
    cap:cfg`participationCap;
    barVol:$[`barVolume in key ctx; ctx`barVolume; 0n];
    maxQty:$[(not null cap) and not null barVol; cap*barVol; 0w];
    filledQty:side*absOrder&maxQty;
    absFilled:abs filledQty;
    / slippage fraction (adverse): bps, optionally * volatility, plus size/ADV impact.
    adv:$[`adv in key ctx; ctx`adv; barVol];
    baseBps:cfg`slippageBps;
    volBps:$[cfg`volScaledSlippage; baseBps*$[`volatility in key ctx; 0f^ctx`volatility; 0f]; baseBps];
    impactBps:$[((cfg`impactCoef)>0f) and (not null adv) and adv>0f; (cfg`impactCoef)*1e4*absFilled%adv; 0f];
    slipFrac:(volBps+impactBps)%1e4;
    proportionalCost:(cfg`proportionalRate)*absFilled;
    slippageCost:absFilled*refPrice*slipFrac;
    fixedCost:$[0f<absFilled; cfg`fixedPerTrade; 0f];
    totalCost:proportionalCost+slippageCost+fixedCost;
    fillPrice:refPrice*1f+side*slipFrac;
    `filledQty`fillPrice`proportionalCost`slippageCost`fixedCost`totalCost!(
        filledQty;fillPrice;proportionalCost;slippageCost;fixedCost;totalCost)
 };

/ ----------------------------------------------------------------------------
/ Execution-realism extensions (v0.77, Research OS R12) - ADDITIVE, opt-in.
/ These are NEW pure functions the replay engine (backtest/replay.q) composes; they
/ do NOT modify .exec.fill. With the replay defaults (.cfg.replay all off/zero) they
/ are no-ops, so a frictionless replay calls the UNCHANGED .exec.fill and reproduces
/ research-mode PnL to tolerance (the canonicals stay byte-identical).
/ ----------------------------------------------------------------------------

/ Participation cap: split a SIGNED order into (filled; remainder) given the as-of bar
/ volume and a participation rate (fraction of volume tradeable in a day). rate<=0 or a
/ null/absent barVolume -> no cap (filled=order, remainder=0). Otherwise the absolute fill
/ is capped at rate*barVolume and the remainder carries the sign. Pure / known-answer testable.
.exec.participationCap:{[order;barVolume;rate]
    noCap:(rate<=0f) or (null barVolume) or (null rate);
    $[noCap;
        `filled`remainder!(order;0f*order);
        [cap:rate*barVolume; filled:(signum order)*(abs order)&cap; `filled`remainder!(filled;order-filled)]]
 };

/ Roll penalty: the extra adverse cost of ROLLING the held exposure on a roll date (R11's roll map
/ flags it). The difference-method return treats a roll as costless (no trade), so the realistic roll
/ cost is charged on the ROLLED EXPOSURE = |position value| * penaltyBps. isRoll 0b / penaltyBps 0f
/ -> 0 (so the frictionless replay is byte-identical). Pure. The replay adds it to the step's cost.
.exec.rollPenalty:{[rolledNotional;isRoll;penaltyBps] $[isRoll; (abs rolledNotional)*penaltyBps%1e4; 0f]};

-1 "execution.q loaded - .exec.* daily fill-and-cost layer ready";
