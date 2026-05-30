/ test_vol_surface_validation.q - surface validation edge cases
\l core/init.q

.test.expectError:{[testName;fn]
    functionResult:@[{x[];`NO_ERROR};fn;{`ERROR}];
    if[functionResult~`NO_ERROR; '"Expected error but got success: ",testName];
    -1 "  PASS ",testName;
 };

/ 1. Empty surface lookup
.test.expectError["empty surface lookup";
    {emptyVolSurface:([] underlying:`symbol$(); expiry:`float$(); strike:`float$(); optionType:`symbol$(); impliedVolatility:`float$());
     .surface.getSurfaceVolatility[emptyVolSurface;100f;1f]}];

/ 2. Surface with negative implied vol
.test.expectError["negative implied vol in surface";
    {badSurface:([] underlying:enlist `AAPL; expiry:enlist 1f; strike:enlist 100f; optionType:enlist `call; impliedVolatility:enlist -0.1);
     .surface.validateVolSurface badSurface}];

/ 3. Surface with zero implied vol
.test.expectError["zero implied vol in surface";
    {badSurface:([] underlying:enlist `AAPL; expiry:enlist 1f; strike:enlist 100f; optionType:enlist `call; impliedVolatility:enlist 0f);
     .surface.validateVolSurface badSurface}];

/ 4. Surface missing required column
.test.expectError["missing impliedVolatility column";
    {badSurface:([] underlying:enlist `AAPL; expiry:enlist 1f; strike:enlist 100f; optionType:enlist `call);
     .surface.validateVolSurface badSurface}];

/ 5. Empty surface validation
.test.expectError["empty surface validation";
    {emptyVolSurface:([] underlying:`symbol$(); expiry:`float$(); strike:`float$(); optionType:`symbol$(); impliedVolatility:`float$());
     .surface.validateVolSurface emptyVolSurface}];

-1 "PASS test_vol_surface_validation";
