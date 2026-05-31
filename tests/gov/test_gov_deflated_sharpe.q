\l core/init.q
/ ============================================================================
/ test_gov_deflated_sharpe.q - PURE-FUNCTION known-answer + direction checks on
/ the PSR / deflated-Sharpe core (.gov.psr / .gov.deflatedSharpe / .gov.phiInv).
/ Synthetic, no HDB. Research OS R3 (gov layer).
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};
apx:{[a;b;tol] tol>abs a-b};

/ --- 1. PSR(SR_hat=0, SR_star=0, large n) = Phi(0) = 0.5 (within the CDF's 7.5e-8) ---
chk[apx[.gov.psr[0f;0f;1000;0f;3f];0.5;1e-6]; "psr(0,0,large n) must be ~0.5"];

/ --- 2. phiInv is the inverse normal: phiInv(0.975) ~ 1.959964 ---
chk[apx[.gov.phiInv 0.975;1.959963985;1e-5]; "phiInv(0.975) must be ~1.96"];
chk[apx[.gov.phiInv 0.5;0f;1e-9]; "phiInv(0.5) must be 0"];

/ --- 3. one hand-computed Phi-argument case ---
/ srHat=0.1, srStar=0, n=101, skew=0, kurt=3 ->
/   num = 0.1*sqrt(100) = 1.0 ; denArg = 1 + ((3-1)/4)*0.1^2 = 1.005 ; arg = 1/sqrt(1.005).
arg:1f%sqrt 1.005;
chk[apx[.gov.psr[0.1;0f;101;0f;3f];.gov.__phi arg;1e-12]; "psr must equal Phi of its hand-computed argument"];
chk[(0.84<.gov.psr[0.1;0f;101;0f;3f]) and 0.845>.gov.psr[0.1;0f;101;0f;3f]; "psr hand case must be ~0.8408"];

/ --- 4. monotonicity: more observations -> higher PSR (for srHat>srStar) ---
chk[.gov.psr[0.1;0f;500;0f;3f] > .gov.psr[0.1;0f;50;0f;3f]; "PSR must rise with n"];

/ --- 5. monotonicity: more TRIALS -> lower DSR (the bar rises with N) ---
chk[.gov.deflatedSharpe[0.15;1000;0f;3f;2;0.0025] > .gov.deflatedSharpe[0.15;1000;0f;3f;50;0.0025]; "DSR must fall as N (trials) rises"];

/ --- 6. monotonicity: positive skew -> higher PSR (favourable, per Bailey & LdP) ---
chk[.gov.psr[0.1;0f;200;0.5;3f] > .gov.psr[0.1;0f;200;0f;3f]; "PSR must rise with positive skew"];

/ --- 7. high per-period Sharpe on a long sample, N=1 -> DSR near 1 ---
chk[.gov.deflatedSharpe[0.15;1000;0f;3f;1;0f] > 0.99; "high SR, long n, N=1 -> DSR ~ 1"];

/ --- 8. the SAME Sharpe with many trials OR a short sample -> DSR well below 0.95 ---
chk[.gov.deflatedSharpe[0.15;1000;0f;3f;50;0.0025] < 0.95; "high SR but N=50 trials -> DSR < 0.95"];
chk[.gov.deflatedSharpe[0.15;20;0f;3f;1;0f] < 0.95; "high SR but only n=20 obs -> DSR < 0.95"];

/ --- 9. degenerate guards: n<=1 -> null PSR ---
chk[null .gov.psr[0.1;0f;1;0f;3f]; "psr with n<=1 must be null"];

-1 "test_gov_deflated_sharpe: PASS";
