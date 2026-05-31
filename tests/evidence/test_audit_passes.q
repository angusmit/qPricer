\l core/init.q
/ ============================================================================
/ test_audit_passes.q - the Evidence Audit passes a CLEAN replay run. R13. SYNTHETIC.
/ A clean synthetic replay run (built by R12's fold over synthetic futures, kept inside the
/ train+validate zone) passes ALL eight evidence checks (overall pass=1b).
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};

mkRows:{[cm;ym;expiryD;fd;ds;px;vol] n:count ds; flip `commodity`contractYM`expiry`firstDate`date`settle`open`high`low`volume!(
    (n#cm);(n#ym);(n#expiryD);(n#fd);ds;px;px;px;px;vol)};
nd:100; dts:2022.01.03+til nd;
e1:dts 55; e2:2022.12.31;
c1dts:dts where dts<=e1; nc1:count c1dts;
c1:mkRows[`SYNE;202203;e1;first dts;c1dts;50f+0.1*`float$til nc1;nc1#5000];
c2:mkRows[`SYNE;202206;e2;first dts;dts;53f+0.08*`float$til nd;nd#5000];
futures:c1,c2;
/ clean run over the first 79 dates (within trainValidate; the boundary is dts[79]).
runDts:79#dts;
rows:flip `stepIndex`stepDate`rawTarget`targetPosition!(til 79; runDts; 79#1f; 79#1f);
foldCfg:`commodity`replayCfg`rollDays`execCfgBase`notional`fromDate`toDate`label!(
    `SYNE;.cfg.replay;5;.exec.__resolve enlist[`txnCostRate]!enlist 0.0005;1f;0Nd;0Nd;`synthLong);
run:.backtest.replay.__fold[rows;foldCfg];

aud:.evidence.audit run;
chk[aud`pass; "a clean replay run passes the evidence audit overall"];
chk[all value aud`checks; "every individual evidence check passes on a clean run"];
chk[(`lookAhead`universe`rollRespected`costsApplied`fillsPresent`bookTies`pnlTies`containment)~key aud`checks; "the audit reports all eight named checks"];
chk[(max run[`steps]`stepDate)<=.evidence.__trainValidateBoundary `SYNE; "the clean run is inside the train+validate zone (containment precondition)"];

-1 "test_audit_passes: PASS";
