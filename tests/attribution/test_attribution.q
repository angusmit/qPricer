\l core/init.q
/ ============================================================================
/ test_attribution.q - PnL explain + bucketed curve risk (.attribution.*). R15. SYNTHETIC.
/ KNOWN-ANSWER decompositions on hand-built runs: a pure PARALLEL move -> level (residual~0); a pure
/ SLOPE move on a deferred position -> slope dominant; a pure ROLL-DOWN -> carry; THE RECONCILIATION
/ (level+slope+curvature+carry+residual == realized total within tol); an UNEXPLAINED move -> high
/ residual fraction; the bucketed risk -> per-tenor delta + roll-down sign.
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};
apx:{[a;b;tol] tol>abs a-b};

/ 6 monthly contracts; build a futures snapshot on `date` with the given per-contract prices.
yms:202406+til 6;
mkSnap:{[date;prices] e:{"D"$(string x div 100),".",(-2#"0",string x mod 100),".20"} each yms; n:6;
    flip `commodity`contractYM`expiry`firstDate`date`settle`open`high`low`volume!(
        (n#`SYN);yms;e;(n#date);(n#date);prices;prices;prices;prices;(n#5000))};
/ a 2-step run holding `heldYM` (position 1) across (prevDate -> curDate); stepPnl[1] = the held
/ contract's realized return (so the decomposition has a real total to reconcile to).
mkRun:{[prevDate;curDate;heldYM;heldRet]
    steps:flip `stepIndex`stepDate`activeContract`position`frontReturn`positionPnl`stepPnl`totalCost`financingCost!(
        (0 1); (prevDate;curDate); (heldYM;heldYM); (1f;1f); (0f;heldRet); (0f;heldRet); (0f;heldRet); (0f;0f); (0f;0f));
    `meta`steps`rollEvents`provenance!((`commodity`fromDate`toDate!(`SYN;prevDate;curDate)); steps; (); ())};
pd:2024.06.03; cd:2024.06.04;
tau:{[date] (((mkSnap[date;6#0f])`expiry)-date)%365};       / per-contract tenor (years) on a date

/ === (1) PURE PARALLEL: flat curve + a uniform +0.5 shift; held front -> level dominates ===
futures:(mkSnap[pd;6#60f]),mkSnap[cd;6#60.5f];
r1:mkRun[pd;cd;first yms;(60.5-60)%60];
a1:.attribution.pnl r1;
chk[a1`reconciled; "(parallel) the decomposition reconciles to the realized total"];
chk[(abs a1`level)>(abs a1`slope)+(abs a1`curvature)+abs a1`carry; "(parallel) level is the dominant component"];
chk[1e-6>a1`residualFraction; "(parallel) the residual fraction is ~0 (fully explained by level)"];

/ === (2) PURE SLOPE (R10 slope shape = k*(tenor-frontTenor)); held a DEFERRED contract -> slope dominates ===
tcd:tau cd; k:6f;
futures:(mkSnap[pd;6#60f]),mkSnap[cd;60f+k*tcd-first tcd];   / curDate = flat + k*(tenor-front)
heldD:yms 4;                                                / a deferred contract (not the front pivot)
heldRet2:((60f+k*(tcd 4)-first tcd)-60f)%60f;               / the deferred contract's realized return
r2:mkRun[pd;cd;heldD;heldRet2];
a2:.attribution.pnl r2;
chk[a2`reconciled; "(slope) the decomposition reconciles"];
chk[(abs a2`slope)>(abs a2`level)|(abs a2`curvature)|abs a2`carry; "(slope) slope is the dominant component on a deferred position"];

/ === (3) PURE ROLL-DOWN: a TENOR-SPACE-STATIC backwardated curve C(tau)=70-8*tau; held front, aged 30
/ days -> the held contract slides UP the static backwardated curve as its tenor shrinks -> carry ===
pd3:2024.06.03; cd3:2024.07.03;                            / 30 days apart -> meaningful roll-down
fr:yms 2;                                                  / a deferred contract that survives both dates
pd3Prices:70f-8f*tau pd3;                                  / each contract priced off the fixed curve C(tau)
cd3Prices:70f-8f*tau cd3;                                  / 30 days later: each contract's tenor shrank
futures:(mkSnap[pd3;pd3Prices]),mkSnap[cd3;cd3Prices];
heldRet3:((cd3Prices 2)-pd3Prices 2)%pd3Prices 2;          / the held contract's realized roll-down return
r3:mkRun[pd3;cd3;fr;heldRet3];
a3:.attribution.pnl r3;
chk[a3`reconciled; "(carry) the decomposition reconciles"];
chk[(abs a3`carry)>(abs a3`level)|(abs a3`slope)|abs a3`curvature; "(carry) carry is the dominant curve-factor component on a static curve"];

/ === (4) UNEXPLAINED: a high-frequency zigzag shift (orthogonal to level/slope/butterfly) -> high residual ===
zig:0.6*(-1 1 -1 1 -1 1f);                                  / alternating, not spanned by the smooth basis
futures:(mkSnap[pd;6#60f]),mkSnap[cd;60f+zig];
heldZ:yms 3;
heldRetZ:((60f+zig 3)-60f)%60f;
r4:mkRun[pd;cd;heldZ;heldRetZ];
a4:.attribution.pnl r4;
chk[a4`reconciled; "(unexplained) the decomposition reconciles"];
chk[0.5<a4`residualFraction; "(unexplained) an out-of-basis move lands in residual with a HIGH residual fraction"];

/ === (5) bucketed risk: a known front position -> delta in the front bucket; roll-down sign matches ===
futures:(mkSnap[pd3;pd3Prices]),mkSnap[cd3;cd3Prices];
rk:.attribution.risk r3;
chk[1f=(rk`bucketDelta)`b0; "(risk) the front position's delta is in the front tenor bucket"];
chk[0f<rk`rollDownExposure; "(risk) roll-down exposure is positive on a backwardated curve (long, +rollYield)"];
chk[1f=rk`calendarSpreadExposure; "(risk) the single-leg calendar-spread exposure equals the position"];

-1 "test_attribution: PASS";
