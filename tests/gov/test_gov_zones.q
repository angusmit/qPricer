\l core/init.q
/ ============================================================================
/ test_gov_zones.q - the train/validate/holdout data zones (.gov.zone.boundaries).
/ Synthetic date vector, no HDB. Research OS R3b (gov layer).
/ Asserts: each zone's [from;to] is right; trainValidate EXCLUDES every holdout date;
/ the zones are contiguous and cover the whole range; holdout is the MOST-RECENT slice.
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};

d:2020.01.01+til 10;                  / 10 contiguous dates
b:.gov.zone.boundaries d;             / default .cfg.gov.zones = 0.6 / 0.2 (holdout 0.2)

/ with 10 dates: nTrain=6, nVal=2, holdout=2.
chk[(b`train)~(2020.01.01;2020.01.06); "train must be the first 6 dates"];
chk[(b`validate)~(2020.01.07;2020.01.08); "validate must be the next 2 dates"];
chk[(b`holdout)~(2020.01.09;2020.01.10); "holdout must be the LAST 2 dates (most recent)"];
chk[(b`trainValidate)~(2020.01.01;2020.01.08); "trainValidate must cover train+validate"];

/ trainValidate EXCLUDES every holdout date (the seal's precondition).
chk[(b[`trainValidate;1]) < b[`holdout;0]; "trainValidate must end strictly before holdout begins"];

/ contiguity + full coverage: train.to+1 = validate.from, validate.to+1 = holdout.from,
/ train.from = min(d), holdout.to = max(d).
chk[(b[`train;1])+1 = b[`validate;0]; "train and validate must be contiguous"];
chk[(b[`validate;1])+1 = b[`holdout;0]; "validate and holdout must be contiguous"];
chk[(b[`train;0]) = min d; "zones must start at the first date"];
chk[(b[`holdout;1]) = max d; "zones must end at the last date (holdout most recent)"];

/ a different fraction split is honoured (40/40/20 over 100 dates -> 40/40/20).
d2:2021.01.01+til 100;
savedZones:.cfg.gov`zones;
.cfg.gov[`zones]:`trainFrac`validateFrac!(0.4;0.4);
b2:.gov.zone.boundaries d2;
chk[(b2`train)~(2021.01.01;2021.01.01+39); "40% train -> first 40 dates"];
chk[(b2`holdout)~(2021.01.01+80;2021.01.01+99); "remaining 20% -> holdout (last 20)"];
.cfg.gov[`zones]:savedZones;          / restore (test independence)

-1 "test_gov_zones: PASS";
