\l core/init.q
/ ============================================================================
/ test_seasonality.q - causal curve/spread seasonality (.season.*). R14. SYNTHETIC.
/ The same-month z-score known-answer; THE KEY INVARIANT - a planted FUTURE same-month observation
/ does NOT change a past z (causal); below min-N the z is null + lowConfidence; the deseasonalised
/ level removes a constant per-calendar-month seasonal factor (known-answer == 1).
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};
apx:{[a;b;tol] tol>abs a-b};

/ --- pure z known-answer ---
chk[apx[0f; .season.__z[3f; 1 2 3 4 5f; 3]; 1e-12]; "z of the mean is 0"];
chk[apx[(5f-3f)%dev 1 2 3 4 5f; .season.__z[5f; 1 2 3 4 5f; 3]; 1e-12]; "z = (x-mean)/dev"];
chk[null .season.__z[5f; 1 2f; 3]; "below min-N the z is null"];

/ --- a synthetic monthly futures with a per-calendar-month seasonal spread + cross-year variation ---
/ monthly contracts 2019.01..2026.12, each trading on the 5th of each month up to its 20th-of-month expiry.
seasLvl:62 60 58 59 61 63 64 62 60 59 61 63f;                  / a NON-linear per-month level (front-M3 spread varies by month)
yms:raze {[y] (y*100)+1+til 12} each 2019+til 8;
mexp:{"D"$(string x div 100),".",(-2#"0",string x mod 100),".20"};
d5:{"D"$(string x div 100),".",(-2#"0",string x mod 100),".05"} each yms;   / the 5th of each contract month
buildF:{[yms;d5;mexp;seasLvl;noiseAmp]
    raze {[ym;d5;mexp;seasLvl;noiseAmp] e:mexp ym; rows:d5 where d5<e; m:ym mod 100;
        lvl:(seasLvl m-1)+noiseAmp*(`float$ym div 100)*(`float$m);   / year x month term: front-M3 spread varies across years
        n:count rows;
        flip `commodity`contractYM`expiry`firstDate`date`settle`open`high`low`volume!(
            (n#`SYN);(n#ym);(n#e);(n#first d5);rows;(n#lvl);(n#lvl);(n#lvl);(n#lvl);(n#5000))}[;d5;mexp;seasLvl;noiseAmp] each yms};

/ dataset A: with cross-year variation -> a finite, meaningful same-month z.
futures:buildF[yms;d5;mexp;seasLvl;0.2];
asOf:2025.01.05;
f1:.season.features[asOf;`SYN];
chk[not null f1`sameMonthZ; "with enough same-month history the z is finite"];
chk[f1[`nSameMonth]>=.cfg.season`minObs; "enough same-calendar-month observations"];

/ THE KEY INVARIANT: a planted FUTURE same-month observation (date>asOf) must NOT change the past z.
future:buildF[enlist 202601;enlist 2026.01.05;mexp;seasLvl;9.9];   / a wildly-different Jan obs in the FUTURE
futures:futures, future;
f2:.season.features[asOf;`SYN];
chk[f1[`sameMonthZ]~f2`sameMonthZ; "a planted FUTURE same-month obs does NOT change a past z (causal/as-of)"];
chk[f1[`seasonalFactor]~f2`seasonalFactor; "the causal seasonal factor is unchanged by future data"];

/ --- min-N: at an EARLY asOf a calendar month has too few observations -> null + lowConfidence ---
futures:buildF[yms;d5;mexp;seasLvl;0.2];
fEarly:.season.features[2019.02.05;`SYN];                  / only 1 January (2019) seen by Feb 2019
chk[fEarly`lowConfidence; "early in the data the same-month history is below min-N (low confidence)"];

/ --- deseasonalised: a CONSTANT per-calendar-month level -> deseasonalised == 1 (known-answer) ---
futures:buildF[yms;d5;mexp;seasLvl;0f];                         / no cross-year variation: each month's front is constant
fc:.season.features[2025.01.05;`SYN];
chk[apx[1f; fc`deseasonalised; 1e-9]; "a constant per-month level deseasonalises to 1 (factor removed)"];

-1 "test_seasonality: PASS";
