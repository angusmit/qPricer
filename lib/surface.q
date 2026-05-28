/ surface.q - simple implied volatility surface: build, lookup, market data
/ Nearest-neighbour lookup in v0.9. No interpolation or arbitrage checks yet.

/ --- Public ---

.surface.buildVolSurface:{[optionChainWithIVs]
    / Filter to successful IV rows only
    okRows:optionChainWithIVs where optionChainWithIVs[`ivStatus]=`OK;
    if[0=count okRows; '"Cannot build vol surface: no successful implied vol rows"];
    ([] underlying:okRows[;`underlying];
        expiry:okRows[;`expiry];
        strike:okRows[;`strike];
        optionType:okRows[;`optionType];
        impliedVolatility:okRows[;`impliedVolatility])
 };

.surface.getSurfaceVolatility:{[volSurface;targetStrike;targetExpiry]
    .surface.getNearestSurfaceVolatility[volSurface;targetStrike;targetExpiry]
 };

.surface.getNearestSurfaceVolatility:{[volSurface;targetStrike;targetExpiry]
    if[0=count volSurface; '"Vol surface is empty"];
    / Find nearest expiry first, then nearest strike within that expiry
    expiryDistances:abs (volSurface`expiry)-targetExpiry;
    nearestExpiryDist:min expiryDistances;
    expiryMatchRows:volSurface where expiryDistances=nearestExpiryDist;
    / Among matched expiry rows, find nearest strike
    strikeDistances:abs (expiryMatchRows`strike)-targetStrike;
    nearestStrikeDist:min strikeDistances;
    matchRows:expiryMatchRows where strikeDistances=nearestStrikeDist;
    / Return first match
    surfaceRow:matchRows 0;
    surfaceRow`impliedVolatility
 };

.surface.createSurfaceMarketData:{[underlying;spot;riskFreeRate;dividendYield;volSurface]
    .surface.validateVolSurface volSurface;
    `underlying`spot`riskFreeRate`dividendYield`marketDataType`volSurface!(
        underlying;spot;riskFreeRate;dividendYield;`volSurface;volSurface)
 };

.surface.validateVolSurface:{[volSurface]
    if[0=count volSurface; '"Vol surface is empty"];
    requiredCols:`expiry`strike`impliedVolatility;
    surfaceCols:cols volSurface;
    missingCols:requiredCols where not requiredCols in surfaceCols;
    if[0<count missingCols; '"Vol surface missing columns: ",", " sv string missingCols];
    if[any (volSurface`impliedVolatility)<=0f; '"Vol surface contains non-positive implied volatility"];
 };
