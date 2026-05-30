\l core/init.q
/ HEADLINE identifiability test: simulate the Schwartz-Smith model with KNOWN
/ params (a specific kappa AND vols), generate a multi-date / multi-tenor futures
/ panel, run the Kalman MLE, and assert the recovered kappa AND volatilities are
/ within tolerance of the known values. This is the proof that the PANEL DYNAMICS
/ identify kappa (and the vols) where v0.51's single curve snapshot could not.
/ ----------------------------------------------------------------------------
/ Tolerances reflect finite-sample estimation error for a ~2-year weekly panel
/ plus the grid resolution; they are far tighter than the snapshot (which cannot
/ identify kappa at all and pins it to a search bound). Deterministic (fixed seed).
trueParams:`kappa`muXi`sigChi`sigXi`correlation`lamChi`lamXi`measSigma!(1.5;0.0;0.25;0.12;0.30;0f;0f;0.008);
taus:0.08 0.17 0.25 0.5 0.75 1.0f;
sim:.commodity.kalman.simulatePanel[trueParams;7;104;taus;2020.01.06;log 60f;42];
panel:sim`panel;
.testutil.assertTrue[104=count distinct panel`obsDate;"panel spans 104 dates"];

est:.commodity.kalman.estimate[panel;()!()];
ep:est`estimatedParams;
.testutil.assertTrue[(est`status)=`OK;"estimation status OK"];

relErr:{[estVal;trueVal] abs (estVal-trueVal)%trueVal};
.testutil.assertTrue[0.20>relErr[ep`kappa;1.5];"kappa identified within 20% (1.5)"];
.testutil.assertTrue[0.20>relErr[ep`sigChi;0.25];"short-factor vol identified within 20% (0.25)"];
.testutil.assertTrue[0.25>relErr[ep`sigXi;0.12];"long-factor vol identified within 25% (0.12)"];
.testutil.assertTrue[0.40>relErr[ep`correlation;0.30];"correlation recovered (looser; weakly identified)"];

/ Identification, not a boundary artifact: kappa sits well inside the [0.05,5] box.
.testutil.assertTrue[((ep`kappa)>0.5) and (ep`kappa)<3.0;"kappa is interior (genuine identification, not a bound)"];
.testutil.assertTrue[not (est`loglik)=0w;"log-likelihood is finite"];
.testutil.assertTrue[104=count est`chi;"filtered chi series one value per date"];
.testutil.assertTrue[104=count est`xi;"filtered xi series one value per date"];

-1 "PASS test_kalman_schwartz_smith_roundtrip: kappa=",(string ep`kappa)," (true 1.5), sigChi=",(string ep`sigChi)," sigXi=",(string ep`sigXi);
