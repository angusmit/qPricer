\l core/init.q
/ ============================================================================
/ test_regime_analogue_distance.q - the analogue distance metric
/ (.regime.analogue.distance). Pure, synthetic fingerprints, no HDB. Research OS R4.
/ Asserts: identical fingerprints -> 0; an axis perturbation gives the expected
/ weighted-Euclidean magnitude + ordering; a categorical mismatch adds the configured penalty.
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};
apx:{[a;b;tol] tol>abs a-b};

k:`curveState`volState`liqState`rollPhase`seasonPhase`slopePct`volPct`volumePct;
a:k!(`backwardation;`high;`deep;`near;`peak;0.9;0.8;0.7);

/ --- identical fingerprints -> 0 ---
chk[0=.regime.analogue.distance[a;a]; "distance to itself must be 0"];

/ --- one percentile axis perturbed by 0.4 (weight 1) -> sqrt(1*0.4^2) = 0.4 ---
b:k!(`backwardation;`high;`deep;`near;`peak;0.5;0.8;0.7);
chk[apx[.regime.analogue.distance[a;b];0.4;1e-12]; "slopePct off by 0.4 (w=1) -> distance 0.4"];

/ --- a nearer perturbation ranks below a farther one (ordering) ---
bNear:k!(`backwardation;`high;`deep;`near;`peak;0.8;0.8;0.7);   / 0.1 off
chk[.regime.analogue.distance[a;bNear] < .regime.analogue.distance[a;b]; "a smaller pct gap -> a smaller distance"];

/ --- one categorical mismatch adds exactly the configured catPenalty (default 0.5) ---
pen:.cfg.regime[`analogue;`catPenalty];
c:k!(`contango;`high;`deep;`near;`peak;0.9;0.8;0.7);            / curveState mismatch only
chk[apx[.regime.analogue.distance[a;c];pen;1e-12]; "one categorical mismatch -> + catPenalty"];

/ --- two categorical mismatches add 2*catPenalty ---
c2:k!(`contango;`low;`deep;`near;`peak;0.9;0.8;0.7);           / curveState + volState mismatch
chk[apx[.regime.analogue.distance[a;c2];2*pen;1e-12]; "two categorical mismatches -> + 2*catPenalty"];

/ --- pct gap AND a categorical mismatch combine (euclid + penalty) ---
mix:k!(`contango;`high;`deep;`near;`peak;0.5;0.8;0.7);          / 0.4 pct gap + 1 cat
chk[apx[.regime.analogue.distance[a;mix];0.4+pen;1e-12]; "pct gap + categorical penalty must add"];

-1 "test_regime_analogue_distance: PASS";
