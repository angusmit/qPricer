\l core/init.q
/ ============================================================================
/ test_workflow_runreplay.q - the R16 end-to-end replay-mode bounded loop (.workflow.runReplay).
/ SYNTHETIC. BOUNDED BY CONSTRUCTION (the R7 discipline): always escalateToHuman, no deploy/allocate/
/ autoTrade field, never tradeable without all gates; the EVIDENCE AUDIT is invoked (a contaminated
/ replay run -> evidenceFailed, the gates never run); the ATTRIBUTION is computed on the run. The
/ existing .workflow.run is unchanged. Gov-state isolation: reset the registry + use UNIQUE hypoIds.
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};
/ GOV-STATE ISOLATION (the R13 caveat): reset the gov tables so registers do not collide with hypos
/ left by the workflow/evidence groups in the full-suite run (the runner values every test in ONE process).
.gov.hypoTbl:.gov.__emptyHypotheses[]; .gov.trialTbl:.gov.__emptyTrials[];
/ CARDS-STATE ISOLATION: the cards tests leave a synthetic `modelCards global; restore the REAL card
/ store (.cfg.cards) so .cards.gateReady finds the seasonalCalSpread card (else runReplay would refuse).
modelCards:.cfg.cards;

/ a synthetic SYNX: monthly contracts over ~9 years with a seasonal-by-delivery-month curve.
mexp:{"D"$(string x div 100),".",(-2#"0",string x mod 100),".20"};
yms:raze {[y] (y*100)+1+til 12} each 2018+til 9;
d5:{"D"$(string x div 100),".",(-2#"0",string x mod 100),".05"} each yms;
seasV:60 60.4 60.9 61.2 61.0 60.5 59.9 59.5 59.6 60.0 60.5 60.9f;
buildC:{[ym;d5;mexp;seasV] e:mexp ym; rows:d5 where d5<e; m:ym mod 100; n:count rows; px:(seasV m-1)+0.02*`float$til n;
    flip `commodity`contractYM`expiry`firstDate`date`settle`open`high`low`volume!(
        (n#`SYNX);(n#ym);(n#e);(n#first d5);rows;px;px;px;px;(n#5000))};
futures:raze buildC[;d5;mexp;seasV] each yms;
cfg:.cfg.strategy.seasonalCalSpread;
/ the run record scoped to train+validate (so the evidence audit's containment check is satisfied).
tv:.gov.zone.range[`SYNX;`trainValidate];
cleanRun:((.template.scs.run `commodity`dateFrom`dateTo!(`SYNX;tv 0;tv 1))`meta)`runRecord;
runner:.template.scs.runner[`SYNX;cfg];
mkProp:{[hid;run] `hypoId`thesis`edgeSource`instruments`claimedRegimes`capName`runRecord`runner!(
    hid;"seasonal cal spread";`mispricing;enlist `SYNX;enlist `backwardation;`seasonalCalSpread;run;.template.scs.runner[`SYNX;.cfg.strategy.seasonalCalSpread])};

/ --- the clean end-to-end loop: a BOUNDED human-escalation packet ---
chk[(.evidence.audit cleanRun)`pass; "(sanity) the trainValidate-scoped run passes the evidence audit"];
pkt:.workflow.runReplay mkProp[`R16_WF_CLEAN;cleanRun];
chk[`escalateToHuman=pkt`decision; "the packet ALWAYS escalates to the human"];
chk[1b=pkt`humanSignOffRequired; "human sign-off is required"];
chk[not any `deploy`allocate`autoTrade in key pkt; "the packet carries NO deploy/allocate/autoTrade field (bounded)"];
chk[-1h=type pkt`anyTradeable; "anyTradeable is a derived boolean fact, not a decision"];
chk[`replay=pkt`mode; "the packet is tagged replay mode"];
chk[all `level`slope`curvature`carry`residual in key pkt`attribution; "the PnL attribution is computed on the run"];

/ --- the EVIDENCE AUDIT is invoked: a contaminated replay run -> evidenceFailed, gates NEVER run ---
s:cleanRun`steps; pd:s`provDateTo; pd[5]:pd[5]+1;          / plant a look-ahead (provDateTo > asOf)
badRun:@[cleanRun;`steps;:;@[s;`provDateTo;:;pd]];
pktBad:.workflow.runReplay mkProp[`R16_WF_BAD;badRun];
chk[`evidenceFailed=(first pktBad`verdicts)`verdict; "a contaminated run yields evidenceFailed (the audit is invoked)"];
chk[0b=pktBad`anyTradeable; "a contaminated run is NOT tradeable (the gates never ran)"];
chk[0b=pktBad`evidencePassed; "evidencePassed is false on the contaminated run"];

/ --- the existing R7 loop is unchanged ---
chk[100h=type .workflow.run; "the existing .workflow.run is intact (unchanged by R16)"];

-1 "test_workflow_runreplay: PASS";
