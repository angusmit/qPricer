\l core/init.q
/ ============================================================================
/ test_factor_rv_replay.q - R17: the replay-mode factor-RV runner (.template.factorRvReplay.*).
/ SYNTHETIC. Asserts: (1) the signal is CAUSAL - a planted FUTURE curve move does NOT change a past
/ signal (the key invariant, inherited from R9's door); (2) the strategy trades the ACTUAL contracts the
/ curve names (the run record PASSES the R13 evidence audit - universe + rollRespected + all 8 checks);
/ (3) R8's factor-stability + residual-stationarity gates still BITE (wired + reject bad inputs); (4) the
/ run through .workflow.runReplay is BOUNDED (escalateToHuman, no deploy field, never tradeable without
/ all gates) and the attribution is computed. Gov-state isolation: reset the registry + UNIQUE hypoIds.
/ R8's template/gates/demo/tests are UNTOUCHED - this reuses them.
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};
/ GOV-STATE ISOLATION (the R13 caveat): the runner values every test in ONE process, so .gov.hypoTbl
/ PERSISTS across groups - reset it + use a UNIQUE hypoId so registers do not collide.
.gov.hypoTbl:.gov.__emptyHypotheses[]; .gov.trialTbl:.gov.__emptyTrials[];
/ CARDS-STATE ISOLATION: the cards tests leave a synthetic `modelCards global; restore the REAL card store
/ (.cfg.cards) so .cards.gateReady finds R8's factorRelativeValue card (else runReplay would refuse).
modelCards:.cfg.cards;

/ ---------------------------------------------------------------------------
/ (1) CAUSALITY: build a levels panel directly and prove a FUTURE move cannot change a PAST signal.
/ A common level factor (a wave) + a contract-specific idiosyncratic residual (the fade-able part).
/ ---------------------------------------------------------------------------
T:90; M:5;
mkCol:{[j;ts] 60f + (0.5*j) + (1.5*sin 0.25*`float$ts) + (0.05*`float$ts) + 0.3*sin (0.7*`float$ts)+j};
levels:flip mkCol[;til T] each til M;                     / T x M
cfgC:.template.factorRvReplay.__cfg[(enlist `lookback)!enlist 20];
pos1:(.template.factorRvReplay.__signal[levels;cfgC])`pos;
/ plant a BIG future curve move at row 75 (>> any past change-row index we assert on).
levels2:levels; levels2[75]:levels[75]+10f;
pos2:(.template.factorRvReplay.__signal[levels2;cfgC])`pos;
tCut:50;                                                 / change-rows < 74 are unaffected; assert at 50
chk[(tCut#pos1)~tCut#pos2; "(1) causal: a planted FUTURE curve move does not change a past signal"];
chk[not (pos1~pos2); "(1) sanity: the future move DID change the later signal (the test is not vacuous)"];

/ ---------------------------------------------------------------------------
/ (3) R8's gates still BITE (reused VERBATIM). residual-stationarity rejects a random walk; factor-
/ stability rejects rotated loadings. (Build before the HDB-backed run so the asserts are self-contained.)
/ ---------------------------------------------------------------------------
rw:sums 1f,(neg 1f),1f,1f,(neg 1f),1f,1f,1f,(neg 1f),1f,1f,(neg 1f),1f,1f,1f,(neg 1f),1f,1f,1f,1f;  / a (drifting) random-walk-ish series
chk[not (.template.rv.stationarity rw)`stationary; "(3) residual-stationarity gate BITES on a non-mean-reverting series"];
h:30;
p1:flip (`float$til h; h#0f; h#0f; h#0f; h#0f);          / variance only in the FRONT column
p2:flip (h#0f; h#0f; h#0f; h#0f; `float$til h);          / variance only in the BACK column
Xbad:p1,p2;                                              / the top loading ROTATES between the two halves
chk[not (.template.factorRv.stability[Xbad;3;.cfg.factor])`stable; "(3) factor-stability gate BITES: rotated loadings -> not stable"];

/ ---------------------------------------------------------------------------
/ Synthetic HDB-backed futures (monthly SYNX) with a common level factor + idiosyncratic residual, so the
/ fade signal actually fires. The HDB query layer reads the global `futures table.
/ ---------------------------------------------------------------------------
mexp:{"D"$(string x div 100),".",(-2#"0",string x mod 100),".20"};
yms:raze {[y] (y*100)+1+til 12} each 2018+til 9;         / 108 monthly contracts
d5:{"D"$(string x div 100),".",(-2#"0",string x mod 100),".05"} each yms;
buildC:{[ym;d5;mexp]
    e:mexp ym; rows:d5 where d5<e; n:count rows; m:ym mod 100;
    di:d5?rows;                                          / global date index of each row
    lvl:5f*sin 0.5*`float$di;                            / common level factor (moves the whole curve)
    idio:0.8*sin (0.9*`float$di)+m;                      / contract-specific idiosyncratic residual
    px:60f + (0.4*m) + lvl + idio;
    flip `commodity`contractYM`expiry`firstDate`date`settle`open`high`low`volume!(
        (n#`SYNX);(n#ym);(n#e);(n#first rows);rows;px;px;px;px;(n#5000))};
futures:raze buildC[;d5;mexp] each yms;
cfg:.template.factorRvReplay.__cfg[(enlist `lookback)!enlist 12];   / small lookback for the short synthetic
tv:.gov.zone.range[`SYNX;`trainValidate];
res:.template.factorRvReplay.run `commodity`dateFrom`dateTo`overrides!(`SYNX;tv 0;tv 1;(enlist `lookback)!enlist 12);
run:(res`meta)`runRecord;

/ (3-wired) the run's validation carries BOTH of R8's gates.
chk[all `factorStability`stationarity in key res`validation; "(3) run validation carries R8's factor-stability + residual-stationarity gates"];
chk[`stable in key res[`validation;`factorStability]; "(3) factor-stability gate is wired into the run"];
chk[`stationary in key res[`validation;`stationarity]; "(3) residual-stationarity gate is wired into the run"];

/ (2) ACTUAL CONTRACTS: the run record passes the R13 evidence audit (universe = live contracts as of the
/ date; rollRespected = forward rolls; + all 8 checks) - the strongest "trades the actual contracts" proof.
chk[0<count run`steps; "(2) the run record has steps"];
aud:.evidence.audit run;
chk[aud`pass; "(2) the run record PASSES the evidence audit (actual contracts + forward rolls + reconciles): ",aud`reason];
acYM:run[`steps]`activeContract;
chk[all acYM in exec contractYM from futures; "(2) every recorded activeContract is a REAL contract id from the curve"];
if[0<count run`rollEvents; chk[all (run[`rollEvents]`toContract)>run[`rollEvents]`fromContract; "(2) every roll is FORWARD in contract"]];

/ (4) BOUNDED end-to-end loop via .workflow.runReplay + the attribution is computed.
prop:`hypoId`thesis`edgeSource`instruments`claimedRegimes`capName`runRecord`runner!(
    `R17_FRV_CLEAN; "factor-RV replay (synthetic)"; `mispricing; enlist `SYNX; `backwardation`contango`flat;
    `factorRelativeValue; run; .template.factorRvReplay.runner[`SYNX;cfg]);
pkt:.workflow.runReplay prop;
chk[`escalateToHuman=pkt`decision; "(4) the packet ALWAYS escalates to the human"];
chk[1b=pkt`humanSignOffRequired; "(4) human sign-off is required"];
chk[not any `deploy`allocate`autoTrade in key pkt; "(4) the packet carries NO deploy/allocate/autoTrade field (bounded)"];
chk[-1h=type pkt`anyTradeable; "(4) anyTradeable is a derived boolean fact, not a decision"];
chk[`replay=pkt`mode; "(4) the packet is tagged replay mode"];
chk[1b=pkt`evidencePassed; "(4) the clean synthetic run passed the evidence audit inside the loop"];
chk[all `level`slope`curvature`carry`residual in key pkt`attribution; "(4) the PnL attribution is computed on the run"];

/ (4-evidence) a CONTAMINATED run -> evidenceFailed, the gates NEVER run (the audit is invoked in the loop).
s:run`steps; pd:s`provDateTo; pd[3]:pd[3]+1;             / plant a look-ahead (provDateTo > asOf)
badRun:@[run;`steps;:;@[s;`provDateTo;:;pd]];
propBad:prop,`hypoId`runRecord!(`R17_FRV_BAD;badRun);
pktBad:.workflow.runReplay propBad;
chk[`evidenceFailed=(first pktBad`verdicts)`verdict; "(4) a contaminated run yields evidenceFailed (the audit is invoked)"];
chk[0b=pktBad`anyTradeable; "(4) a contaminated run is NOT tradeable (the gates never ran)"];

/ R8 is UNTOUCHED: its template function is intact and distinct from the replay runner.
chk[100h=type .template.factorRv.run; "R8's factorRelativeValue template is intact (unchanged by R17)"];
chk[not .template.factorRv.run~.template.factorRvReplay.run; "R17's replay runner is a DISTINCT function (R8 not rebased)"];

-1 "test_factor_rv_replay: PASS";
