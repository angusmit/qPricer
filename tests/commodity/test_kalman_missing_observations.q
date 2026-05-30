\l core/init.q
/ The Kalman filter handles a varying observation set per date (missing tenors):
/ build a panel, drop some (date,tenor) cells, and assert the filter still runs,
/ produces a finite log-likelihood, and yields one filtered state per date.
trueParams:`kappa`muXi`sigChi`sigXi`correlation`lamChi`lamXi`measSigma!(1.2;0.0;0.20;0.10;0.20;0f;0f;0.01);
taus:0.08 0.25 0.5 0.75 1.0f;
sim:.commodity.kalman.simulatePanel[trueParams;14;16;taus;2020.01.06;log 55f;7];
panel:sim`panel;
fullDates:count distinct panel`obsDate;

/ Drop the 0.5y tenor on the 3rd date and the 1.0y tenor on the 8th date.
dropDates:(distinct panel`obsDate)1 7;
sparsePanel:select from panel where not (
    ((obsDate=dropDates 0) and tau=0.5) or ((obsDate=dropDates 1) and tau=1.0));
.testutil.assertTrue[(count sparsePanel)<count panel;"some observations removed"];
.testutil.assertTrue[fullDates=count distinct sparsePanel`obsDate;"all dates still present"];

filtered:.commodity.kalman.filter[sparsePanel;trueParams];
.testutil.assertTrue[(filtered`status)=`OK;"filter runs on the sparse panel"];
.testutil.assertTrue[not (filtered`loglik)=0w;"log-likelihood finite with missing observations"];
.testutil.assertTrue[not (filtered`loglik)=-0w;"log-likelihood not negative-infinite"];
.testutil.assertTrue[fullDates=count filtered`chi;"one filtered chi per date despite missing tenors"];

-1 "PASS test_kalman_missing_observations: loglik=",string filtered`loglik;
