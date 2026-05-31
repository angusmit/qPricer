\l core/init.q
/ ============================================================================
/ test_curve_snapshot.q - immutable daily curve snapshots (.curve.snapshot*). Research OS R10.
/ Synthetic snapshot rows (no HDB). A snapshot written then read back is identical; re-writing
/ the same (commodity, date) is a no-op-if-matching (immutability); a DIFFERING re-write errors
/ (catches a non-deterministic curve); the snapshots table is a single flat (unpartitioned) table.
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};
caught:{[f;arg] @[f;arg;{[e] 1b}]~1b};

/ a synthetic built-curve object -> a snapshot row.
mkCurve:{[px] ([] tenor:0.1*1+til count px; price:px; contractYM:202003+til count px; expiry:2020.03.20+30*til count px)};
mkObj:{[asOf;comm;effDate;px] cv:mkCurve px;
    `asOf`commodity`effectiveDate`curve`spreads`features!(asOf;comm;effDate;cv;.curve.__spreads cv;.curve.__features[cv;.cfg.curve])};

.curve.snapshots:.curve.__emptySnap[];

/ --- write a snapshot, read it back ---
row:.curve.__snapRow mkObj[2020.02.03;`SYN;2020.02.03;60 59 58 57 56 55f];
.curve.__snapWrite row;
back:.curve.snapshotAt[2020.02.03;`SYN];
chk[1=count back; "one snapshot row written"];
chk[(first back)~row; "the snapshot reads back identical to what was written"];
chk[((first back)`classification)~`backwardation; "the snapshot carries the derived classification"];

/ --- re-snapshotting the SAME (commodity,date) with a MATCHING curve is a no-op (immutable) ---
.curve.__snapWrite .curve.__snapRow mkObj[2020.02.03;`SYN;2020.02.03;60 59 58 57 56 55f];
chk[1=count .curve.snapshotAt[2020.02.03;`SYN]; "a matching re-snapshot is a no-op (no duplicate, no overwrite)"];

/ --- a DIFFERING re-write for the same (commodity,date) ERRORS (catches non-determinism) ---
chk[caught[.curve.__snapWrite; .curve.__snapRow mkObj[2020.02.03;`SYN;2020.02.03;55 56 57 58 59 60f]];
    "a differing snapshot for the same (commodity,date) must be REJECTED (immutability)"];
chk[((first .curve.snapshotAt[2020.02.03;`SYN])`classification)~`backwardation; "the rejected write did not corrupt the existing snapshot"];

/ --- a second date appends; history is ordered and unpartitioned (one flat table) ---
.curve.__snapWrite .curve.__snapRow mkObj[2020.03.02;`SYN;2020.03.02;55 56 57 58 59 60f];
h:.curve.snapshotHistory[`SYN];
chk[2=count h; "a new (commodity,date) appends a second snapshot"];
chk[(h`date)~asc h`date; "snapshot history is date-ordered"];
chk[98h=type .curve.snapshots; "curveSnapshots is a single flat (unpartitioned) table"];

-1 "test_curve_snapshot: PASS";
