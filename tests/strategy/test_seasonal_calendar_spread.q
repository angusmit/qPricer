\l core/init.q
/ ============================================================================
/ test_seasonal_calendar_spread.q - the R16 capstone template (.template.scs.*). SYNTHETIC.
/ The signal is correct (z>+thr -> SHORT the spread, z<-thr -> LONG, exit on z-cross-of-zero / maxHold);
/ the causal seasonal z uses only same-calendar-month obs BEFORE t; the STATIONARITY GATE bites
/ (rejects a random-walk spread, passes a mean-reverting one). PRE-REGISTERED params, not tuned.
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};
apx:{[a;b;tol] tol>abs a-b};

/ --- the position state machine (pure, known-answer) ---
z1:(0n;0n;2.0;0.5;-0.3;-2.0;1.0);                          / null,null,wide,hold,revert,narrow,revert
pos1:.template.scs.__positions[z1;1.5;20];
chk[pos1~(0;0;-1;-1;0;1;0); "signal: z>+entry SHORTs the spread, z<-entry LONGs, exit on z-cross-of-zero"];
/ maxHold: a persistently-wide z exits at maxHold then can re-enter.
pos2:.template.scs.__positions[5#2.0;1.5;3];
chk[pos2~(-1;-1;-1;0;-1); "signal: a persistently-wide spread exits at maxHold, then re-enters"];
chk[0=first .template.scs.__positions[enlist 0n;1.5;20]; "a null z (below min-N) gives no signal (stays flat)"];

/ --- the causal seasonal z (same-calendar-month, strictly before t; min-N gated) ---
dts:2020.01.15+0 31 62 365 396 427;                       / Jan,Feb,Mar 2020; Jan,Feb,Mar 2021
spread:1.0 2.0 3.0 5.0 6.0 7.0f;
z:.template.scs.__seasonalZ[spread;dts;1];                / minObs 1
chk[null z 0; "the first observation has no prior same-month history -> null"];
chk[null z 1; "Feb 2020: no prior February -> null"];
/ Jan 2021 (idx 3): prior same-month (January) = {1.0}; one obs, dev 0 -> null (need variation).
chk[null z 3; "Jan 2021 vs a single prior January (zero dispersion) -> null"];

/ --- the stationarity gate BITES (reuse .template.rv.stationarity on the spread) ---
mrSpread:{[n] s:0f; r:(); do[n; s:(0.6*s)+0.3*((count r) mod 5)-2; r,:s]; r}[60];  / OU-ish, AR1~0.6
rwSpread:10f+0.5*`float$til 60;                            / a linear trend: AR(1)~1, non-stationary
chk[(.template.rv.stationarity mrSpread)`stationary; "stationarity gate PASSES a mean-reverting spread"];
chk[not (.template.rv.stationarity rwSpread)`stationary; "stationarity gate REJECTS a random-walk spread (fading noise)"];

-1 "test_seasonal_calendar_spread: PASS";
