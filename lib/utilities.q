/ utilities.q — validation helpers, interpolation, error measurement

/ --- Validation ---

.utilities.requireKeys:{[dict;k;name]
    missing:k where not k in key dict;
    if[count missing; '"Missing required field in ",name,": ",string first missing];
 };

.utilities.assertPositive:{[val;name]
    if[not val > 0; 'name," must be positive"];
 };

.utilities.assertNonNegative:{[val;name]
    if[val < 0; 'name," must be non-negative"];
 };

.utilities.assertSupportedValue:{[val;supported;name]
    if[not val in supported;
        '"Unsupported ",name,": ",string[val],". Supported values: ","," sv string supported];
 };

.utilities.assertBoolean:{[val;name]
    if[not val in 01b; 'name," must be a boolean (0b or 1b)"];
 };

.utilities.assertPositiveInteger:{[val;name]
    if[not (val > 0) and val = floor val; 'name," must be a positive integer"];
 };

/ --- Interpolation ---

.utilities.findNearestIndex:{[xVal;xGrid]
    (-1 + count xGrid) & 0 | -1 + first where xGrid >= xVal
 };

.utilities.linearInterpolate:{[xVal;xGrid;yGrid]
    if[(xVal < first xGrid) or xVal > last xGrid;
        '"Interpolation failed: value ",string[xVal]," is outside grid range [",string[first xGrid],", ",string[last xGrid],"]"];
    if[xVal >= last xGrid; :last yGrid];
    i:.utilities.findNearestIndex[xVal;xGrid];
    if[xGrid[i] = xVal; :yGrid i];
    w:(xVal - xGrid i) % xGrid[i+1] - xGrid i;
    yGrid[i] + w * yGrid[i+1] - yGrid i
 };

/ --- Error measurement ---

.utilities.relativeError:{[actual;expected]
    $[expected = 0f; abs actual; abs (actual - expected) % expected]
 };
