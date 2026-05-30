/ seasonality.q - deterministic, opt-in seasonal overlay for commodity curves (v0.51)
/ ----------------------------------------------------------------------------
/ A pure deterministic seasonal adjustment applied ADDITIVELY to the log-forward
/ of a commodity curve: F_seasonal(t) = F(t) * exp(factor(t)). The factor is a
/ smooth annual cycle by default, or an explicit 12-month vector if supplied.
/ ----------------------------------------------------------------------------
/ Most relevant to gas / power (winter demand peaks, summer troughs); crude is
/ only mildly seasonal. This overlay is OPT-IN on the curve adapters
/ (.strategy.path.fromFuturesCurve / fromCorrelatedCurves) via an optional
/ seasonCfg key - when absent, those adapters are byte-identical to before. It
/ will be exercised on natural-gas (NG) curves once that data is added.
/ ----------------------------------------------------------------------------
/ seasonCfg keys:
/   amplitude   - peak log-forward adjustment (e.g. 0.10 = +/-10% at the extremes)
/   phaseYears  - calendar-year fraction at which the cycle peaks (0 = year start)
/   monthlyFactors (optional) - explicit 12-vector of log adjustments by month;
/                 when present it OVERRIDES the smooth cosine cycle.
/ factor[timeYears; seasonCfg] is vectorised over a timeYears vector. timeYears is
/ the delivery time in years (its fractional part selects the season).

.commodity.seasonality.__twoPi:6.28318530717958647692;

.commodity.seasonality.factor:{[timeYears;seasonCfg]
    if[`monthlyFactors in key seasonCfg;
        monthlyFactors:`float$seasonCfg`monthlyFactors;
        if[not 12=count monthlyFactors; '"seasonality: monthlyFactors must be a 12-vector"];
        fracYear:timeYears-floor timeYears;
        monthIdx:floor 12f*fracYear;
        :monthlyFactors monthIdx];
    amplitude:`float$seasonCfg`amplitude;
    phaseYears:`float$seasonCfg`phaseYears;
    amplitude*cos .commodity.seasonality.__twoPi*timeYears-phaseYears
 };

/ Fit a 12-vector of monthly seasonal factors from observed (timeYears, value)
/ pairs by averaging the values within each calendar month (the fractional part
/ of timeYears selects the month). The recovered vector is exactly the
/ monthlyFactors a seasonal overlay would use. Months with no observation -> null.
.commodity.seasonality.fitMonthlyFactors:{[timeYears;values]
    if[not (count timeYears)=count values; '"fitMonthlyFactors: timeYears/values length mismatch"];
    months:floor 12f*timeYears-floor timeYears;
    {[months;values;m] monthValues:values where months=m; $[0=count monthValues;0n;avg monthValues]}[months;`float$values;] each til 12
 };
