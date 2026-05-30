\l core/init.q
/ ============================================================================
/ kalman_schwartz_smith.q - real-data example (NOT a test).
/ ----------------------------------------------------------------------------
/ Part B on the REAL WTI panel: build a futures panel across 2020-2021 (non-
/ positive prices excluded), run the Kalman-filter MLE, and print the estimated
/ Schwartz-Smith params with an IDENTIFIED kappa (from the curve DYNAMICS, which
/ a single snapshot could not), the filtered chi/xi factor series, and the
/ implied spot level. CSVs are user-supplied and gitignored.
/ ============================================================================
crudeDir:"data/barchart/CRUDE";
longTable:.parser.crude.loadAll crudeDir;
allDates:asc distinct longTable`date;
inRange:allDates where (allDates>=2020.01.01) and allDates<=2021.12.31;
asofDates:inRange where 0=(til count inRange) mod 21;
hist:.parser.crude.curveHistory[longTable;asofDates];
panel:.commodity.kalman.panelFromCurveHistory hist;
-1 "Panel: ",(string count distinct panel`obsDate)," dates, ",(string count panel)," observations.";
-1 "Running Kalman MLE (coordinate descent) ...";

est:.commodity.kalman.estimate[panel;()!()];
ep:est`estimatedParams;
-1 "";
-1 "Kalman-MLE Schwartz-Smith parameters (kappa IDENTIFIED from the panel dynamics):";
-1 "  kappa (mean reversion) = ",string ep`kappa;
-1 "  sigChi (short vol)     = ",string ep`sigChi;
-1 "  sigXi  (long vol)      = ",string ep`sigXi;
-1 "  rho                    = ",string ep`correlation;
-1 "  muXi   (equil drift)   = ",string ep`muXi;
-1 "  measSigma              = ",string ep`measSigma;
-1 "  logLik                 = ",string est`loglik;
-1 "";
-1 "Filtered factor series (chi = short-term deviation, xi = equilibrium level):";
show ([] obsDate:est`dates; chi:est`chi; xi:est`xi);
-1 "";
-1 "Implied spot exp(chi+xi) at last date: ",string exp (last est`chi)+last est`xi;
-1 "Versus v0.51 (single snapshot) which pinned kappa to a bound, the panel here";
-1 "identifies kappa; feeding this kappa into Part A gives a cleaner cy series.";
