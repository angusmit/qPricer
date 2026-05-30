/ test_vol_surface_nearest_lookup.q - nearest-neighbour surface lookup
\l core/init.q

volSurface:([]
    underlying:`AAPL`AAPL`AAPL`AAPL;
    expiry:0.5 0.5 1.0 1.0;
    strike:90 100 100 110f;
    optionType:`call`call`call`call;
    impliedVolatility:0.22 0.20 0.21 0.23);

/ 1. Exact lookup: strike=100, expiry=0.5
exactVol:.surface.getSurfaceVolatility[volSurface;100f;0.5];
if[(abs exactVol-0.20)>0.001; '"FAIL: exact lookup expected 0.20, got ",string exactVol];

/ 2. Near strike: strike=102, expiry=0.5 -> nearest strike is 100
nearStrikeVol:.surface.getSurfaceVolatility[volSurface;102f;0.5];
if[(abs nearStrikeVol-0.20)>0.001; '"FAIL: near-strike lookup expected 0.20, got ",string nearStrikeVol];

/ 3. Near expiry then strike: strike=108, expiry=0.9
/ Nearest expiry to 0.9 is 1.0. Within expiry=1.0, nearest strike to 108 is 110.
nearExpiryVol:.surface.getSurfaceVolatility[volSurface;108f;0.9];
if[(abs nearExpiryVol-0.23)>0.001; '"FAIL: near-expiry lookup expected 0.23, got ",string nearExpiryVol];

/ 4. All returned vols positive
if[not exactVol>0f; '"FAIL: exact vol not positive"];
if[not nearStrikeVol>0f; '"FAIL: near-strike vol not positive"];
if[not nearExpiryVol>0f; '"FAIL: near-expiry vol not positive"];

-1 "PASS test_vol_surface_nearest_lookup: exact=",string[exactVol],", nearStrike=",string[nearStrikeVol],", nearExpiry=",string nearExpiryVol;
