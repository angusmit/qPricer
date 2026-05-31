\l core/init.q
/ ============================================================================
/ test_audit_bites.q - the Evidence Audit BITES: one deliberately-contaminated run PER CHECK, each
/ caught with the right check failing. R13. SYNTHETIC. An audit that always passes is worse than
/ useless - so every check must fail the run it should (the same discipline as R2/R5/R6).
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};
/ amend one column of a run record's steps table (functional - NO qSQL on a local table).
setCol:{[run;col;v] @[run;`steps;:;@[run`steps;col;:;v]]};
failed:{[run;c] not (.evidence.audit run)[`checks][c]};   / 1b iff check c failed on this run

mkRows:{[cm;ym;expiryD;fd;ds;px;vol] n:count ds; flip `commodity`contractYM`expiry`firstDate`date`settle`open`high`low`volume!(
    (n#cm);(n#ym);(n#expiryD);(n#fd);ds;px;px;px;px;vol)};
nd:100; dts:2022.01.03+til nd;
e1:dts 55; e2:2022.12.31;
c1dts:dts where dts<=e1; nc1:count c1dts;
c1:mkRows[`SYNE;202203;e1;first dts;c1dts;50f+0.1*`float$til nc1;nc1#5000];
c2:mkRows[`SYNE;202206;e2;first dts;dts;53f+0.08*`float$til nd;nd#5000];
futures:c1,c2;
foldCfg:`commodity`replayCfg`rollDays`execCfgBase`notional`fromDate`toDate`label!(
    `SYNE;.cfg.replay;5;.exec.__resolve enlist[`txnCostRate]!enlist 0.0005;1f;0Nd;0Nd;`synthLong);
runDts:79#dts;
rows:flip `stepIndex`stepDate`rawTarget`targetPosition!(til 79; runDts; 79#1f; 79#1f);
run:.backtest.replay.__fold[rows;foldCfg];
chk[(.evidence.audit run)`pass; "(sanity) the base run is clean before contamination"];
s:run`steps;

/ 1. LOOK-AHEAD: one step accesses data dated AFTER its asOf (provDateTo > asOf).
pd:s`provDateTo; pd[5]:pd[5]+1;
chk[failed[setCol[run;`provDateTo;pd];`lookAhead]; "look-ahead check catches a provDateTo > asOf"];

/ 2. UNIVERSE: a fill on a contract not live as of the date (a fabricated contractYM).
ac:s`activeContract; ac[0]:999999;
chk[failed[setCol[run;`activeContract;ac];`universe]; "universe check catches a fill on a non-live contract"];

/ 3. ROLL-RESPECTED: a rollEvents row that disagrees with the active-contract transitions.
bogus:([] date:enlist runDts 70; commodity:enlist `SYNE; fromContract:enlist 202206j; toContract:enlist 202203j; reason:enlist `bogus);
chk[failed[@[run;`rollEvents;:;(run`rollEvents),bogus];`rollRespected]; "roll-respected check catches a fabricated roll event"];

/ 4. COSTS-APPLIED: a non-frictionless run with the transaction cost zeroed out.
chk[failed[setCol[run;`proportionalCost;(count s)#0f];`costsApplied]; "costs-applied check catches missing costs on a run with fills"];

/ 5. FILLS-PRESENT: a position change with no corresponding fill.
pos:s`position; pos[10]:pos[10]+5f;
chk[failed[setCol[run;`position;pos];`fillsPresent]; "fills-present check catches a position change with no fill"];

/ 6. BOOK-TIES: the position book offset from the cumulative fills.
chk[failed[setCol[run;`position;10f+s`position];`bookTies]; "book-ties check catches a book that does not equal cumulative fills"];

/ 7. PnL-TIES: a stepPnl that does not reconcile to positionPnl - costs - financing.
sp:s`stepPnl; sp[5]:sp[5]+1f;
chk[failed[setCol[run;`stepPnl;sp];`pnlTies]; "PnL-ties check catches PnL that does not reconcile"];

/ 8. CONTAINMENT: a run whose dates span the holdout zone.
rowsFull:flip `stepIndex`stepDate`rawTarget`targetPosition!(til nd; dts; nd#1f; nd#1f);
runFull:.backtest.replay.__fold[rowsFull;foldCfg];
chk[failed[runFull;`containment]; "containment check catches a run that spans the holdout"];

-1 "test_audit_bites: PASS";
