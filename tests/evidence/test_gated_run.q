\l core/init.q
/ ============================================================================
/ test_gated_run.q - the fail-safe precondition wrapper .evidence.gatedRun. R13. SYNTHETIC.
/ On a FAILED audit it returns `evidenceFailed and the gates NEVER run (the runner is not invoked);
/ on a PASSED audit it delegates to the gates (.gov.runFull, which invokes the runner). The existing
/ .workflow.run is unchanged. Mirrors .cards.gatedRun (audit -> on-pass gate / on-fail reject).
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};

mkRows:{[cm;ym;expiryD;fd;ds;px;vol] n:count ds; flip `commodity`contractYM`expiry`firstDate`date`settle`open`high`low`volume!(
    (n#cm);(n#ym);(n#expiryD);(n#fd);ds;px;px;px;px;vol)};
nd:120; dts:2022.01.03+til nd;
e1:dts 55; e2:2023.06.30;
c1dts:dts where dts<=e1; nc1:count c1dts;
c1:mkRows[`SYNG;202203;e1;first dts;c1dts;50f+0.1*`float$til nc1;nc1#5000];
c2:mkRows[`SYNG;202206;e2;first dts;dts;53f+0.08*`float$til nd;nd#5000];
futures:c1,c2;
foldCfg:`commodity`replayCfg`rollDays`execCfgBase`notional`fromDate`toDate`label!(
    `SYNG;.cfg.replay;5;.exec.__resolve enlist[`txnCostRate]!enlist 0.0005;1f;0Nd;0Nd;`synthLong);
runDts:90#dts;                                               / within trainValidate (boundary dts[95])
rows:flip `stepIndex`stepDate`rawTarget`targetPosition!(til 90; runDts; 90#1f; 90#1f);
cleanRun:.backtest.replay.__fold[rows;foldCfg];
chk[(.evidence.audit cleanRun)`pass; "(sanity) the clean run passes the audit"];

/ a contaminated run (planted look-ahead) the audit will reject.
s:cleanRun`steps; pd:s`provDateTo; pd[5]:pd[5]+1;
badRun:@[cleanRun;`steps;:;@[s;`provDateTo;:;pd]];
chk[not (.evidence.audit badRun)`pass; "(sanity) the contaminated run fails the audit"];

/ A counting runner so we can observe whether the gates ran. We use an UNREGISTERED hypo id so the
/ test stays ISOLATED from shared gov state (no .gov.register, which in the full-suite value-context
/ would upsert into a hypoTbl already populated by other groups) - .gov.runFull raises for an
/ unregistered hypo, and that it RAISES (rather than returning the evidenceFailed short-circuit) is
/ exactly the proof that the wrapper took the delegation branch.
runnerCalls:0;
runner:{[from;to] runnerCalls+:1; sub:dts where dts within (from;to); ([] date:sub; pnl:(count sub)#0.001)};

/ FAIL path: gatedRun on the contaminated run -> evidenceFailed, and the gates NEVER run.
runnerCalls:0;
vFail:.evidence.gatedRun[`H_UNREG_EV;badRun;runner;`curveState];
chk[`evidenceFailed=(first vFail)`verdict; "a failed audit yields the evidenceFailed verdict"];
chk[0b=(first vFail)`tradeable; "the evidenceFailed verdict is not tradeable"];
chk[0=runnerCalls; "on a failed audit the gates NEVER run (the runner was never invoked)"];

/ PASS path: a passing audit -> gatedRun delegates to .gov.runFull. The hypo is deliberately
/ unregistered, so .gov.runFull raises ("no hypothesis"); the wrapper RAISING (caught here as
/ `delegated) rather than returning the evidenceFailed short-circuit proves it took the delegation
/ branch past the audit gate. (The full successful delegation is shown in the demo, in a fresh process.)
delegated:@[{.evidence.gatedRun[`H_UNREG_EV;x;runner;`curveState];`ran}; cleanRun; {`delegated}];
chk[delegated~`delegated; "on a passed audit the wrapper delegates to the gates (not the evidenceFailed short-circuit)"];

/ the existing R7 workflow loop is intact (unchanged by R13).
chk[100h=type .workflow.run; "the existing .workflow.run is intact (unchanged)"];

-1 "test_gated_run: PASS";
