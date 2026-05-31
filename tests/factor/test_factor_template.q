\l core/init.q
/ ============================================================================
/ test_factor_template.q - the factorRelativeValue template. Research OS R8.
/ The template-specific gates BITE: FACTOR STABILITY passes a stable structure and fails a
/ regime break (the top factor rotates level->slope); RESIDUAL STATIONARITY passes a mean-
/ reverting deviation and fails a random walk. The run composes the RV engine faithfully on
/ the factor deviation. Synthetic, no HDB.
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};
apx:{[a;b;tol] tol>abs a-b};

cfg:.cfg.factor;
Llevel:.factor.__norm 1 1 1 1f;
Lslope:.factor.__norm -3 -1 1 3f;
nObs:160; h:nObs div 2;
s1:0.03*sin 0.3*til nObs;
s2:0.02*cos 0.4*til nObs;

/ --- FACTOR STABILITY gate ---
/ stable: both halves are level-dominated (same top factor) -> high cosine -> PASS.
Xstable:(s1 *\: Llevel) + s2 *\: Lslope;
chk[(.template.factorRv.stability[Xstable;cfg`k;cfg])`stable; "a stable factor structure PASSES the stability gate"];
/ regime break: first half pure level, second half pure slope -> top factor rotates -> FAIL.
w1:nObs#1f; w1[h+til nObs-h]:0f;
w2:nObs#0f; w2[h+til nObs-h]:1f;
Xbreak:((s1*w1) *\: Llevel) + (s2*w2) *\: Lslope;
sb:.template.factorRv.stability[Xbreak;cfg`k;cfg];
chk[not sb`stable; "a regime break (top factor rotates level->slope) FAILS the stability gate"];
chk[(sb`cos)<cfg`minLoadingStability; "the regime-break loading cosine is below the bar"];

/ --- RESIDUAL STATIONARITY gate (reuses the RV AR(1) test on the traded deviation) ---
arDev:{[n] e:n#0.2 -0.3 0.1 -0.15 0.25 -0.2; r:0f; o:(); i:0; while[i<n; r:(0.5*r)+e i; o,:r; i+:1]; o}[150];
rwDev:sums 150#0.1 -0.1 0.2 -0.2 0.15 -0.05;
chk[(.template.rv.stationarity arDev)`stationary; "a mean-reverting residual PASSES the stationarity gate"];
chk[not (.template.rv.stationarity rwDev)`stationary; "a random-walk residual FAILS the stationarity gate"];

/ --- end-to-end run on the stable synthetic panel -> well-formed output ---
out:.template.factorRv.run[`commodity`panel`dates`overrides!(
    `SYN; Xstable; 2020.01.01+til nObs; `lookback`entryZ`txnCostRate!(20;1.0;0f))];
chk[(`pnl`validation`meta)~key out; "run returns pnl/validation/meta"];
chk[nObs=count out`pnl; "one PnL row per change-panel date"];
chk[(out[`meta;`shape])~`factorRelativeValue; "meta records the factorRelativeValue shape"];
chk[(out[`meta;`edgeSource])~`structural; "the template claims a STRUCTURAL edge"];
chk[`factorStability in key out`validation; "validation carries the factor-stability gate"];
chk[`stationarity in key out`validation; "validation carries the residual-stationarity gate"];
chk[(cfg`k)=count out[`meta;`explainedVar]; "meta reports explained variance per PC (k of them)"];

/ --- faithful composition: the PnL is the RV engine applied to the factor deviation ---
d:.factor.pca[Xstable;cfg`k;cfg];
devSrs:sums (d`residuals)[;cfg`residualMaturity];
expected:select date, pnl from .template.rv.__signalPnl[devSrs; 2020.01.01+til nObs; .cfg.templates.rv,`lookback`entryZ`txnCostRate!(20;1.0;0f)];
chk[(out`pnl)~expected; "the template composes the RV mean-reversion engine on the factor deviation (faithful)"];

-1 "test_factor_template: PASS";
