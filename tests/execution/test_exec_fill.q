\l core/init.q
/ Known-answer unit tests for the .exec.fill daily fill-and-cost model: cost
/ components computed independently, participation cap partial fill, adverse
/ slippage direction, and the byte-identical legacy default (proportional only).

/ buy turnover 10 at ref 50, slippage 20 bps, propRate 5e-4, fixed 0.1, no cap.
cfgBuy:.cfg.exec,`proportionalRate`slippageBps`fixedPerTrade!(0.0005;20f;0.1);
rb:.exec.fill[10f;`refPrice`currentPos!(50f;0f);cfgBuy];
.testutil.assertNear[rb`filledQty;10f;1e-12;"full fill when no cap"];
.testutil.assertNear[rb`proportionalCost;0.0005*10f;1e-12;"proportionalCost = rate*|filled|"];
.testutil.assertNear[rb`slippageCost;10f*50f*20f%1e4;1e-12;"slippageCost = |filled|*ref*bps/1e4"];
.testutil.assertNear[rb`fixedCost;0.1;1e-12;"fixedCost charged when trading"];
.testutil.assertNear[rb`totalCost;(rb`proportionalCost)+(rb`slippageCost)+rb`fixedCost;1e-12;"components add up"];
.testutil.assertTrue[(rb`fillPrice)>50f;"buy fill price is adverse (above ref)"];

/ sell turnover -4 at ref 50, slippage 10 bps: adverse fill price is BELOW ref.
rs:.exec.fill[-4f;`refPrice`currentPos!(50f;10f);.cfg.exec,(enlist `slippageBps)!enlist 10f];
.testutil.assertNear[rs`filledQty;-4f;1e-12;"sell filled qty signed"];
.testutil.assertNear[rs`slippageCost;4f*50f*10f%1e4;1e-12;"sell slippageCost on |filled|"];
.testutil.assertTrue[(rs`fillPrice)<50f;"sell fill price is adverse (below ref)"];

/ participation cap: order 10, bar volume 30, cap 0.1 -> maxQty 3 -> partial fill 3.
rc:.exec.fill[10f;`refPrice`barVolume`currentPos!(50f;30f;0f);.cfg.exec,(enlist `participationCap)!enlist 0.1];
.testutil.assertNear[rc`filledQty;3f;1e-12;"participation cap limits fill to cap*volume"];
/ order 2 under the same cap (maxQty 3) fills fully.
rc2:.exec.fill[2f;`refPrice`barVolume`currentPos!(50f;30f;0f);.cfg.exec,(enlist `participationCap)!enlist 0.1];
.testutil.assertNear[rc2`filledQty;2f;1e-12;"order below cap fills fully"];

/ no trade -> no fixed cost.
rz:.exec.fill[0f;`refPrice`currentPos!(50f;5f);cfgBuy];
.testutil.assertNear[rz`fixedCost;0f;1e-12;"no fixed cost when nothing trades"];
.testutil.assertNear[rz`totalCost;0f;1e-12;"no cost when nothing trades"];

/ LEGACY DEFAULT: zero slippage / fixed, no cap -> totalCost == proportionalRate*|order|.
rd:.exec.fill[8f;`refPrice`currentPos!(50f;0f);.cfg.exec];
.testutil.assertNear[rd`slippageCost;0f;1e-12;"default has zero slippage"];
.testutil.assertNear[rd`totalCost;(.cfg.exec`proportionalRate)*8f;1e-12;"default totalCost == legacy proportional cost"];

-1 "PASS test_exec_fill: known-answer fill + cost components + cap + legacy default";
