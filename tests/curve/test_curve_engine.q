\l core/init.q
/ ============================================================================
/ test_curve_engine.q - the curve engine (.curve.*). Research OS R10. Synthetic, no HDB.
/ Classification/slope/roll-yield/curvature on known curves; the three shock operators; the
/ regime AGREEMENT cross-check (curve classification == regime/'s on the same curve); determinism.
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};
apx:{[a;b;tol] tol>abs a-b};
cfg:.cfg.curve;
mk:{[px] ([] tenor:0.1*1+til count px; price:px; contractYM:202003+til count px; expiry:2020.03.20+30*til count px)};

back:mk 60 59 58 57 56 55f;        / decreasing -> backwardation
cont:mk 55 56 57 58 59 60f;        / increasing -> contango
flat:mk 60 60.01 60.005 60.02 60.01 60.015f;

/ --- classification + slope sign + roll-yield sign ---
fb:.curve.__features[back;cfg]; fc:.curve.__features[cont;cfg]; ff:.curve.__features[flat;cfg];
chk[(fb`classification)~`backwardation; "decreasing curve -> backwardation"];
chk[(fc`classification)~`contango; "increasing curve -> contango"];
chk[(ff`classification)~`flat; "near-flat curve -> flat"];
chk[(fb`slope)<0f; "backwardation slope is negative ((deferred-front)/front)"];
chk[(fc`slope)>0f; "contango slope is positive"];
chk[(fb`rollYield)>0f; "roll yield is positive in backwardation (ln(front/deferred)>0)"];
chk[(fc`rollYield)<0f; "roll yield is negative in contango"];

/ --- curvature known-answer: near - 2*mid + far on 10 12 11 = 10 - 24 + 11 = -3 ---
chk[apx[-3f; (.curve.__features[mk 10 12 11f;cfg])`curvature; 1e-12]; "butterfly curvature = near - 2*mid + far"];

/ --- the three shock operators on tenors 0, 0.5, 1 ---
sc:([] tenor:0 0.5 1f; price:10 10 10f; contractYM:1 2 3; expiry:3#2020.06.01);
chk[((.curve.shock.parallel[sc;1f])`price)~11 11 11f; "parallel shock adds a constant to every tenor"];
chk[((.curve.shock.slope[sc;1f])`price)~10 10.5 11f; "slope shock adds a tenor-proportional amount (front pivot)"];
chk[((.curve.shock.butterfly[sc;1f])`price)~10 11 10f; "butterfly shock bows the belly, leaving the ends fixed"];

/ --- AGREEMENT with regime/: same classification on the same front/deferred (no rebase yet) ---
panel:{[front;deferred] ([] date:2020.01.10 2020.01.13; frontSettle:2#front; frontVolume:2#100f;
    frontExpiry:2#2020.03.20; frontTenorDays:2#70i; deferredSettle:2#deferred)};
agree:{[c]
    di:.cfg.regime`deferredIdx;
    cv:.curve.__features[c;.cfg.curve];
    rp:.regime.__labelPanel[panel[first c`price; (c`price) di]; .cfg.regime];
    (cv`classification)~last rp`curveState};
chk[agree back; "curve classification AGREES with regime/ on a backwardated curve"];
chk[agree cont; "curve classification AGREES with regime/ on a contango curve"];
chk[agree flat; "curve classification AGREES with regime/ on a flat curve"];
chk[(.cfg.curve`flatThreshold)~.cfg.regime`flatSlopeThreshold; "the curve flat threshold is SOURCED FROM regime (cannot drift)"];
chk[(.cfg.curve`deferredIdx)~.cfg.regime`deferredIdx; "the curve deferred index is sourced from regime"];

/ --- determinism: two builds of the same input match byte-for-byte ---
chk[(.curve.__features[back;cfg])~.curve.__features[back;cfg]; "features are deterministic"];

-1 "test_curve_engine: PASS";
