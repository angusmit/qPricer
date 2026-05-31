\l core/init.q
/ ============================================================================
/ test_replay_execution.q - the execution-realism extensions on the replay engine. R12. SYNTHETIC.
/ The participation cap limits a fill to a fraction of the as-of bar volume (remainder carries/drops
/ per config); the roll penalty charges the rolled exposure on a roll date; and with realism OFF the
/ fills equal the UNCHANGED .exec.fill output (frictionless equivalence -> byte-identical canonicals).
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};
apx:{[a;b;tol] tol>abs a-b};

/ --- unit: .exec.participationCap (pure, known-answer) ---
cp:.exec.participationCap[2f;100f;0.001];                     / cap = 0.001*100 = 0.1
chk[apx[0.1; cp`filled; 1e-12]; "participationCap caps the fill at rate*barVolume"];
chk[apx[1.9; cp`remainder; 1e-12]; "participationCap leaves the unfilled remainder"];
cn:.exec.participationCap[-2f;100f;0.001];
chk[apx[-0.1; cn`filled; 1e-12]; "participationCap preserves the order sign"];
chk[(2f~(.exec.participationCap[2f;100f;0f])`filled) and 0f~(.exec.participationCap[2f;100f;0f])`remainder; "rate<=0 -> no cap (full fill, zero remainder)"];
chk[2f~(.exec.participationCap[2f;0n;0.5])`filled; "a null bar volume -> no cap"];

/ --- unit: .exec.rollPenalty (pure, known-answer) ---
chk[apx[0.5; .exec.rollPenalty[100f;1b;50f]; 1e-12]; "rollPenalty charges |notional|*bps on a roll date"];
chk[apx[0f; .exec.rollPenalty[100f;0b;50f]; 1e-12]; "rollPenalty is zero off a roll date"];
chk[apx[0f; .exec.rollPenalty[100f;1b;0f]; 1e-12]; "rollPenalty zero when penaltyBps is zero (frictionless)"];
chk[apx[0.5; .exec.rollPenalty[-100f;1b;50f]; 1e-12]; "rollPenalty uses |notional| (sign-independent)"];

/ --- a synthetic single-contract futures (no roll) for the cap + frictionless-equivalence checks ---
mkRows:{[cm;ym;expiryD;fd;ds;px;vol] n:count ds; flip `commodity`contractYM`expiry`firstDate`date`settle`open`high`low`volume!(
    (n#cm); (n#ym); (n#expiryD); (n#fd); ds; px; px; px; px; vol)};
nd:16; dts:2022.01.03+til nd;
px:50f+(0.2*`float$til nd)+0.6*`float$(til nd) mod 3;
futures:mkRows[`SYN;202212;2022.12.20;first dts;dts;px;nd#100];   / volume 100, far expiry -> no roll
notional:1f; txnRate:0.0005;
/ alternating long/short target -> a large order every day.
rows:flip `stepIndex`stepDate`rawTarget`targetPosition!(til nd; dts; nd#1f; notional*(-1 1) (til nd) mod 2);
baseFold:{[rows;rc;cm] `commodity`replayCfg`rollDays`execCfgBase`notional`fromDate`toDate`label!(
    cm; rc; 5; .exec.__resolve enlist[`txnCostRate]!enlist 0.0005; 1f; 0Nd; 0Nd; cm)};

/ frictionless: every fill is full (filled==order), no slippage/roll cost (== unchanged .exec.fill).
rrF:.backtest.replay.__fold[rows; baseFold[rows;.cfg.replay;`SYN]];
sF:rrF`steps;
chk[all (sF`filledQty)=sF`order; "realism OFF: every order fills fully (== the unchanged .exec.fill)"];
chk[0f=max abs sF`slippageCost; "realism OFF: zero slippage"];
chk[0f=max abs sF`rollCost; "realism OFF: zero roll cost"];

/ participation cap ON: cap = 0.001*100 = 0.1; the big alternating orders are clipped to +/-0.1.
capCfg:.cfg.replay,`participationRate`remainderPolicy!(0.001;`drop);
rrC:.backtest.replay.__fold[rows; baseFold[rows;capCfg;`SYN]];
sC:rrC`steps;
chk[(max abs sC`filledQty) <= 0.1+1e-12; "participation cap ON: no fill exceeds rate*barVolume"];
chk[any (abs sC`filledQty) < abs sC`order; "participation cap ON: the cap actually binds (filled < order)"];

/ --- roll penalty integration: a 2-contract constant-long replay, penalty bites once at the roll ---
e1:2022.03.18; e2:2022.06.18; n2:80; d2:2022.01.03+til n2;
c1d:d2 where d2<=e1; nc1:count c1d;
c1:mkRows[`SYN2;202203;e1;first d2;c1d;50f+0.1*`float$til nc1;nc1#9000];
c2:mkRows[`SYN2;202206;e2;first d2;d2;53f+0.08*`float$til n2;n2#9000];
futures:c1,c2;
rows2:flip `stepIndex`stepDate`rawTarget`targetPosition!(til n2; d2; n2#1f; n2#1f);
rrNoPen:.backtest.replay.__fold[rows2; baseFold[rows2;.cfg.replay;`SYN2]];
penCfg:.cfg.replay,enlist[`rollPenaltyBps]!enlist 50f;
rrPen:.backtest.replay.__fold[rows2; baseFold[rows2;penCfg;`SYN2]];
sPen:rrPen`steps;
rollRow:first select from sPen where isRoll;
expectedRollCost:(abs rollRow[`position]*rollRow`refPrice)*50f%1e4;
/ NOTE: position here is post-fill; on a pure roll the order is 0 so pre==post; the penalty uses the held exposure.
chk[0<rollRow`rollCost; "roll penalty ON: a positive roll cost is charged on the roll date"];
chk[1=count select from sPen where isRoll; "exactly one roll date in the 2-contract series"];
chk[0f=max abs (sPen`rollCost) - ?[sPen`isRoll; sPen`rollCost; 0f]; "roll cost is charged ONLY on roll dates"];
chk[((rrNoPen`meta)`totalPnl) > (rrPen`meta)`totalPnl; "the roll penalty degrades net PnL (cost drag)"];
chk[apx[((rrNoPen`meta)`totalPnl)-(rrPen`meta)`totalPnl; sum sPen`rollCost; 1e-12]; "net degradation equals the total roll cost charged"];

-1 "test_replay_execution: PASS";
