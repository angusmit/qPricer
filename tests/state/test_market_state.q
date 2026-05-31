\l core/init.q
/ ============================================================================
/ test_market_state.q - the Market State object (.state.build). Research OS R9.
/ Assembles the expected dict (curve / spreads / features-placeholder / universe / refs /
/ provenance) from a synthetic futures table; the universe excludes a contract expired as of
/ asOf; two as-of dates produce states where the earlier is a strict time-prefix of the later.
/ No real HDB.
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};

/ synthetic futures: three live contracts trading on the effective date, PLUS an already-
/ expired contract (202001, expiry 2020.01.20) that traded earlier (so it is in the as-of
/ slice but must NOT be in the tradable universe as of 2020.02.03).
futures:([] commodity:`SYN`SYN`SYN`SYN`SYN`SYN`SYN;
    contractYM:202001 202003 202006 202009 202003 202006 202009;
    expiry:2020.01.20 2020.03.20 2020.06.20 2020.09.20 2020.03.20 2020.06.20 2020.09.20;
    firstDate:7#2019.11.01;
    date:2020.01.10 2020.02.03 2020.02.03 2020.02.03 2020.01.10 2020.01.10 2020.01.10;
    settle:48 52 51 50 51.5 50.5 49.5f; open:7#0f; high:7#0f; low:7#0f; volume:7#100);

asOf:2020.02.03;
st:.state.build[asOf;`SYN];

/ --- the dict has the expected contract fields ---
chk[(`asOf`commodity`effectiveDate`curve`spreads`features`universe`refs`provenance)~key st; "state has the expected fields"];
chk[(st`asOf)=asOf; "state records the asOf"];
chk[(st`commodity)~`SYN; "state records the commodity"];
chk[(st`effectiveDate)=2020.02.03; "effective date is the latest trading date <= asOf"];
chk[98h=type st`curve; "curve is a table (the as-of curve)"];
chk[0<count st`curve; "the as-of curve is non-empty"];
chk[98h=type st`spreads; "spreads is a table (adjacent calendar spreads)"];
chk[((count st`curve)-1)=count st`spreads; "one calendar spread per adjacent maturity pair"];
chk[(0=count st`features); "features is a declared (empty) placeholder for R10/R14"];

/ --- the universe excludes the expired contract (202001) and includes the live ones ---
chk[not 202001 in st`universe; "an expired contract (expiry<=asOf) is EXCLUDED from the universe"];
chk[202003 in st`universe; "a live contract is in the universe"];
chk[all `costModel`rollRule`commodity in key st`refs; "refs carries the cost-model + roll-rule handles"];

/ --- the universe-live invariant holds (the as-of curve carries no expired contract) ---
chk[.state.invariant.universeLive st; "the as-of live curve carries no contract expired by its effective date"];

/ --- two as-of dates: the earlier is a strict TIME-PREFIX of the later ---
stEarly:.state.build[2020.01.10;`SYN];
chk[(stEarly`effectiveDate) < st`effectiveDate; "an earlier asOf -> an earlier (strictly prior) effective date"];
chk[((stEarly`provenance)`dateTo) <= ((st`provenance)`dateTo); "the earlier state's data is a time-prefix of the later"];
chk[((stEarly`provenance)`nRows) <= ((st`provenance)`nRows); "the earlier state sees no more rows than the later"];

-1 "test_market_state: PASS";
