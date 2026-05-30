\l core/init.q
/ .strategy.__portfolioValue is a pure helper: PV = cash + legMarkSum + hedgePosition*spot.
/ Strategy-agnostic, used as the independent-revaluation oracle in accounting tests.

pv1:.strategy.__portfolioValue[100f;-50f;0.5;200f];
.testutil.assertNear[pv1;150f;1e-12;"PV = 100 + (-50) + 0.5*200 = 150"];

pv2:.strategy.__portfolioValue[0f;0f;0f;100f];
.testutil.assertNear[pv2;0f;1e-12;"empty book: PV = 0"];

pv3:.strategy.__portfolioValue[-100f;120f;-1f;110f];
.testutil.assertNear[pv3;-90f;1e-12;"borrowed cash with short hedge: PV = -100 + 120 + (-1)*110 = -90"];

pv4:.strategy.__portfolioValue[50f;30f;0f;75f];
.testutil.assertNear[pv4;80f;1e-12;"no hedge: PV = cash + legMarkSum"];

pv5a:.strategy.__portfolioValue[100f;-20f;0.3;101f];
pv5b:.strategy.__portfolioValue[100f;-20.5;0.3;102f];
deltaPV:pv5b-pv5a;
expectedFromComponents:(-0.5f)+0.3f*(102f-101f);
.testutil.assertNear[deltaPV;expectedFromComponents;1e-12;"delta-PV equals legMark change + hedge*spot change when cash unchanged"];

-1 "PASS test_strategy_portfolio_value_helper: pv1=",string[pv1],", deltaPV=",string[deltaPV];
