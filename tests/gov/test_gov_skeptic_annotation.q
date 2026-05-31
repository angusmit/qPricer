\l core/init.q
/ ============================================================================
/ test_gov_skeptic_annotation.q - the regime risk-memory SKEPTIC (v0.71 wiring).
/ The gov verdict gains an informational `riskMemory annotation, but the skeptic changes
/ NO gate pass/fail and NO previously-existing verdict field. Synthetic (no HDB), so the
/ annotation is empty (graceful) - the point here is that the gate logic is byte-identical.
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};

.gov.hypoTbl:.gov.__emptyHypotheses[];
.gov.trialTbl:.gov.__emptyTrials[];
.gov.register `hypoId`thesis`edgeSource`instruments`claimedRegimes`status!(
    `hSk;"synthetic";`riskPremium;`SYN;enlist `backwardation;`research);

rets:0.0008+0.0004*sin 0.3*til 200;
dates:2020.01.01+til 200;
pnl:([] date:dates; pnl:rets);
labs:([] date:dates; curveState:200#`backwardation);

v:.gov.run[`hSk;pnl;labs;`curveState];
chk[1=count v; "one regime bucket -> one verdict"];
row:first v;

/ --- the verdict GAINS a riskMemory annotation ---
chk[`riskMemory in cols v; "the verdict must carry a `riskMemory annotation field"];
chk[(row`riskMemory)~""; "synthetic (no HDB / no library) -> empty annotation, graceful"];

/ --- and the gate logic is BYTE-IDENTICAL: every pre-existing field equals a direct
/ .gov.evaluate (the skeptic added nothing to the cascade). N=1 trial was logged by the run. ---
fam:.gov.trials `hSk;
annDays:.cfg.gov`annualizationDays;
varSR:$[1<count fam; var (fam`netSharpe)%sqrt annDays; 0f];
expect:.gov.evaluate `hypo`rets`bucket`axis`nTrials`varSR`cfg!(
    .gov.hypo `hSk; rets; `backwardation; `curveState; count fam; varSR; .cfg.gov);
chk[(row`verdict)~expect`verdict; "the skeptic changes NO verdict (pass/fail unchanged)"];
chk[(row`failedGate)~expect`failedGate; "failedGate unchanged"];
chk[(row`tradeable)~expect`tradeable; "tradeable unchanged"];
chk[(row`postHoc)~expect`postHoc; "postHoc unchanged"];
chk[(row`netSharpe)~expect`netSharpe; "netSharpe unchanged"];
chk[(row`dsr)~expect`dsr; "dsr unchanged"];

-1 "test_gov_skeptic_annotation: PASS";
