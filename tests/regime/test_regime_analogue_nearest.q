\l core/init.q
/ ============================================================================
/ test_regime_analogue_nearest.q - the analogue ranking (.regime.analogue.nearest /
/ .regime.analogue.rank) on a SYNTHETIC in-memory episode set. No CSV/HDB read.
/ Research OS R4. Asserts: the closest episode ranks first and a clearly-distant one
/ last; the n cut is honoured; the regimeEpisodes schema round-trips (label/driversKey/
/ riskMemory carried through alongside the computed distance).
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};

ek:`episodeId`commodity`dateFrom`dateTo`label`driversKey`riskMemory`curveState`volState`liqState`rollPhase`seasonPhase`slopePct`volPct`volumePct;
mkEp:{[ek;id;cs;sl] ek!(id;`CRUDE;2020.01.01;2020.02.01;id;`drv;"risk memory for ",string id;cs;`high;`deep;`near;`peak;sl;0.8;0.7)};

/ build a synthetic episode set (set the global the library reads).
regimeEpisodes:(upsert/)[.regime.library.__empty[];(
    mkEp[ek;`near_ep;`backwardation;0.88];      / closest to the query below
    mkEp[ek;`mid_ep;`backwardation;0.55];        / a moderate pct gap
    mkEp[ek;`far_ep;`contango;0.05])];           / categorical mismatch + big pct gap

k:`curveState`volState`liqState`rollPhase`seasonPhase`slopePct`volPct`volumePct;
query:k!(`backwardation;`high;`deep;`near;`peak;0.9;0.8;0.7);

res:.regime.analogue.nearest[query;3];
chk[3=count res; "nearest must return all 3 episodes when n=3"];
chk[(first res`episodeId)~`near_ep; "the closest episode must rank FIRST"];
chk[(last res`episodeId)~`far_ep; "the clearly-distant episode must rank LAST"];

/ distances are sorted ascending.
chk[(res`distance)~asc res`distance; "results must be sorted by ascending distance"];

/ the n cut is honoured.
chk[1=count .regime.analogue.nearest[query;1]; "n=1 must return exactly the single nearest"];
chk[(first (.regime.analogue.nearest[query;1])`episodeId)~`near_ep; "n=1 returns the nearest episode"];

/ the schema round-trips: identity + drivers + risk memory carried alongside distance.
chk[all `episodeId`commodity`label`driversKey`riskMemory`distance in cols res; "nearest must carry episode identity + risk memory + distance"];
chk[(first res`riskMemory)~"risk memory for near_ep"; "the joined risk-memory summary must follow its episode"];

/ empty library -> empty result (no crash).
regimeEpisodes:.regime.library.__empty[];
chk[0=count .regime.analogue.nearest[query;3]; "an empty library yields no analogues"];

-1 "test_regime_analogue_nearest: PASS";
