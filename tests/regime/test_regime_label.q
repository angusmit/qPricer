\l core/init.q
/ Market State Engine known-answers on a synthetic futures long table (no CSV, no HDB).
/ A front contract (near expiry) sits below... above a deferred contract -> backwardation;
/ an injected volatile stretch -> volState `high; a thin-volume window -> liqState `thin;
/ percentiles in [0,1]; a hand-checked rolling percentile.
n:60;
dts:2020.01.02+til n;
/ front settle: calm 60-ish, with a volatile stretch in the middle (large alternating moves).
frontS:60f+0.02*til n;
frontS[30+til 8]:frontS[30+til 8]+5f*(-1 1 -1 1 -1 1 -1 1f);    / high-vol burst
defS:frontS-2f;                                                  / deferred 2 BELOW front -> backwardation
frontVol:n#1000f; frontVol[15+til 6]:40f;                        / thin-liquidity window
colNames:`contractYM`expiry`firstDate`date`settle`open`high`low`volume;
front:flip colNames!(n#202012; n#2020.12.20; n#2019.01.01; dts; frontS; frontS; frontS; frontS; frontVol);
defc:flip colNames!(n#202106; n#2021.06.20; n#2019.01.01; dts; defS; defS; defS; defS; n#2000f);
lt:front,defc;
cfg:.regime.defaultConfig[];

panel:.regime.__panelFromLong[lt;cfg];
.testutil.assertTrue[n=count panel;"one panel row per date"];
.testutil.assertTrue[(panel`frontSettle)~frontS;"front = smallest-tenor contract settle"];
.testutil.assertTrue[(panel`deferredSettle)~defS;"deferred = next-tenor contract settle"];

labs:.regime.__labelPanel[panel;cfg];
.testutil.assertTrue[(`date`curveState`volState`liqState`rollPhase`seasonPhase`slopePct`volPct`volumePct)~cols labs;"regime table schema"];
.testutil.assertTrue[all (labs`curveState)=`backwardation;"front>deferred everywhere -> backwardation"];
.testutil.assertTrue[`high in labs`volState;"injected volatile stretch -> volState high somewhere"];
.testutil.assertTrue[`thin in labs`liqState;"thin-volume window -> liqState thin somewhere"];
.testutil.assertTrue[all (labs[`slopePct]>=0f)&labs[`slopePct]<=1f;"slopePct in [0,1]"];
.testutil.assertTrue[all (labs[`volPct]>=0f)&labs[`volPct]<=1f;"volPct in [0,1]"];
.testutil.assertTrue[all (labs[`volumePct]>=0f)&labs[`volumePct]<=1f;"volumePct in [0,1]"];

/ contango known-answer: flip the curve (deferred ABOVE front) -> contango.
panelC:.regime.__panelFromLong[(flip colNames!(n#202012;n#2020.12.20;n#2019.01.01;dts;frontS;frontS;frontS;frontS;frontVol)),flip colNames!(n#202106;n#2021.06.20;n#2019.01.01;dts;frontS+3f;frontS+3f;frontS+3f;frontS+3f;n#2000f);cfg];
.testutil.assertTrue[all (.regime.__labelPanel[panelC;cfg]`curveState)=`contango;"deferred>front everywhere -> contango"];

/ hand-checked rolling percentile (fraction of trailing window <= current; floor-free).
.testutil.assertNear[(last .regime.__rollPct[3;1 2 3 4 5f]);1f;1e-12;"rollPct of an increasing series: last == 1"];
.testutil.assertNear[(.regime.__rollPct[3;5 4 3 2 1f])4;1%3;1e-12;"rollPct w=3 at last of [5 4 3 2 1]: window (3 2 1), 1<=1 -> 1/3"];

-1 "PASS test_regime_label: curve/vol/liquidity states + percentiles known-answer correct";
