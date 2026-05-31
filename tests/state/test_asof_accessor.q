\l core/init.q
/ ============================================================================
/ test_asof_accessor.q - the single as-of door (.state.asof). Research OS R9.
/ On a SYNTHETIC futures table: across several asOf dates the accessor returns NO row with
/ date>asOf; a deliberately-planted FUTURE row is excluded; the provenance stamp matches
/ exactly the rows returned; and the revision-aware path (a synthetic table WITH validDate)
/ takes the latest validDate<=asOf. No real HDB.
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};

/ synthetic futures (the last row is a deliberately-planted FUTURE observation).
futures:([] commodity:`SYN`SYN`SYN`SYN`SYN;
    contractYM:202003 202006 202003 202006 202009;
    expiry:2020.03.20 2020.06.20 2020.03.20 2020.06.20 2020.09.20;
    firstDate:5#2019.12.01;
    date:2020.01.10 2020.01.10 2020.02.15 2020.02.15 2020.07.01;   / 2020.07.01 = planted future
    settle:50 49 52 51 55f; open:5#0f; high:5#0f; low:5#0f; volume:5#100);

/ (i) across several asOf dates: NEVER a row with date>asOf.
{[asOf] a:.state.asof[asOf;`SYN];
    chk[.state.invariant.asofRespected[asOf;a`data]; "accessor returned a row with date>asOf for asOf=",string asOf];
    chk[.state.invariant.provenanceConsistent a; "provenance inconsistent for asOf=",string asOf]
 } each 2020.01.05 2020.01.31 2020.03.01 2020.06.01 2020.12.31;

/ (iii) the planted FUTURE row (2020.07.01) is excluded for an asOf before it.
a:.state.asof[2020.06.01;`SYN];
chk[not 2020.07.01 in (a`data)`date; "the planted future row must be excluded"];
chk[4=count a`data; "exactly the 4 rows with date<=2020.06.01 are returned"];

/ (iv) the provenance stamp matches exactly what was returned.
chk[((a`provenance)`nRows)=count a`data; "provenance nRows must equal the rows returned"];
chk[((a`provenance)`dateTo)=max (a`data)`date; "provenance dateTo must equal the latest returned date"];
chk[((a`provenance)`asOf)=2020.06.01; "provenance records the asOf"];
chk[(asc (a`provenance)`contracts)~asc 202003 202006; "provenance lists exactly the returned contracts"];

/ the append-only provenance log grew by these calls.
chk[0<count .state.provLog; "the provenance log is appended to"];

/ revision-aware path: a synthetic table WITH validDate -> latest revision validDate<=asOf.
futures:([] commodity:`R`R`R; date:3#2020.01.05;
    validDate:2020.01.06 2020.01.20 2020.02.10; px:10 11 12f);
r:.state.asof[2020.01.31;`R];
chk[1=count r`data; "the revision path returns one row per observation date (latest revision)"];
chk[11f=first (r`data)`px; "the revision path takes the latest validDate<=asOf (11, validDate 2020.01.20)"];
chk[not 2020.02.10 in (r`data)`validDate; "a revision published AFTER asOf is excluded"];

-1 "test_asof_accessor: PASS";
