\l core/init.q
/ ============================================================================
/ test_replay_loop.q - the event-driven replay fold (.backtest.replay.__fold). R12. SYNTHETIC.
/ The fold is DETERMINISTIC (two runs identical); a frictionless replay of a synthetic constant-long
/ futures strategy reproduces an INDEPENDENT vectorised computation to tolerance; the run record's PnL
/ TIES to positions x price-changes + costs (the R13 reconciliation); the book ties to the fills; and
/ NO step used data with date>asOf (the as-of invariant, via the per-step provenance).
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};
apx:{[a;b;tol] tol>abs a-b};

/ --- synthetic futures: 2 monthly contracts, c1 rolls to c2 ~5d before its expiry ---
nd:80; dts:2022.01.03+til nd;
e1:2022.03.18; e2:2022.06.18;
c1dts:dts where dts<=e1; nc1:count c1dts;
c1px:50f+(0.1*`float$til nc1)+0.5*`float$(til nc1) mod 3;     / trending + wiggle -> varied returns
c2px:53f+(0.08*`float$til nd)+0.4*`float$(til nd) mod 4;
mkRows:{[cm;ym;expiryD;fd;ds;px;vol] n:count ds; flip `commodity`contractYM`expiry`firstDate`date`settle`open`high`low`volume!(
    (n#cm); (n#ym); (n#expiryD); (n#fd); ds; px; px; px; px; vol)};
c1:mkRows[`SYN;202203;e1;first dts;c1dts;c1px;nc1#5000];
c2:mkRows[`SYN;202206;e2;first dts;dts;c2px;nd#5000];
futures:c1,c2;

/ --- a constant-long futures strategy: target = notional every day ---
notional:1f; txnRate:0.0005;
rows:flip `stepIndex`stepDate`rawTarget`targetPosition!(til nd; dts; nd#1f; nd#notional);
foldCfg:`commodity`replayCfg`rollDays`execCfgBase`notional`fromDate`toDate`label!(
    `SYN; .cfg.replay; 5; .exec.__resolve enlist[`txnCostRate]!enlist txnRate; notional; 0Nd; 0Nd; `synthLong);
rr:.backtest.replay.__fold[rows;foldCfg];
steps:rr`steps;
chk[nd=count steps; "the fold emits one step per trading date"];

/ --- determinism: a second identical fold reproduces the run record exactly ---
rr2:.backtest.replay.__fold[rows;foldCfg];
chk[(rr`steps)~rr2`steps; "the fold is deterministic (two runs are byte-identical)"];

/ --- the active contract + roll: exactly one c1 -> c2 roll, where d2e(c1) crosses 5 ---
chk[1=count rr`rollEvents; "exactly one roll event (c1 -> c2)"];
chk[202203=first (rr`rollEvents)`fromContract; "rolls FROM the front contract c1"];
chk[202206=first (rr`rollEvents)`toContract; "rolls TO the next contract c2"];

/ --- faithfulness to an INDEPENDENT vectorised computation ---
/ the active contract per date (days_before_expiry, decided as of the PRIOR date - the lagged roll the
/ research path uses) + its own return.
priceK:`date`contractYM xkey select date,contractYM,settle from futures;
refDates:(first dts),-1_dts;                                 / refDate[i]=dates[i-1]; refDate[0]=dates[0]
activeRef:?[(e1-refDates)>5; 202203; 202206];
chk[(steps`activeContract)~activeRef; "the fold's active contract matches the independent (prior-date) roll rule"];
pNow:priceK[([] date:1_dts; contractYM:1_activeRef)]`settle;
pPrev:priceK[([] date:-1_dts; contractYM:1_activeRef)]`settle;
refRet:(pNow%pPrev)-1f;
/ constant-long: total = notional*sum(active returns from t>=1) - the t0 entry cost (frictionless).
expectedTotal:(notional*sum refRet)-txnRate*notional;
chk[apx[expectedTotal; (rr`meta)`totalPnl; 1e-10]; "frictionless replay reproduces the vectorised total PnL to tolerance"];
chk[1e-12>max abs refRet-1_steps`frontReturn; "the re-derived as-of returns match the reference returns elementwise"];

/ --- reconciliation (the R13 checks) ---
recon:(steps`stepPnl)-((steps`positionPnl)-((steps`totalCost)+steps`financingCost));
chk[1e-12>max abs recon; "PnL reconciles: stepPnl == positionPnl - totalCost - financingCost"];
posByMove:(0f^prev steps`position)*steps`frontReturn;
chk[1e-12>max abs (steps`positionPnl)-posByMove; "positionPnl == prior position x price-move (ties to positions+moves)"];
chk[1e-12>max abs (steps`position)-sums steps`filledQty; "the book ties to the fills (position == cumulative filled qty)"];

/ --- the as-of invariant: no step used data with date>asOf ---
prov:rr`provenance;
chk[not any (prov`provDateTo)>prov`asOf; "NO step used data with date>asOf (point-in-time safe)"];

-1 "test_replay_loop: PASS";
