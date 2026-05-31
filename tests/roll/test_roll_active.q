\l core/init.q
/ ============================================================================
/ test_roll_active.q - as-of-only roll selection (.roll.active). Research OS R11. Synthetic.
/ days_before_expiry selects the right active contract + rolls exactly rollDays before expiry,
/ with the multi-day window blending old->new; THE KEY INVARIANT: a planted FUTURE volume spike
/ does NOT change a past roll decision (as-of-only); volume_switch rolls when as-of volume
/ crosses; oi_switch flags unavailable.
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};
apx:{[a;b;tol] tol>abs a-b};
caught:{[f;arg] @[f;arg;{[e] 1b}]~1b};

/ synthetic 2-contract futures: c1 (202003, expiry 2020.03.20), c2 (202006, expiry 2020.06.20).
d:2020.01.01+til 130;
c1:([] commodity:`SYN; contractYM:202003; expiry:2020.03.20; firstDate:2020.01.01; date:d; open:130#0f; high:130#0f; low:130#0f; settle:50f+0.05*til 130; volume:130#100);
c2:([] commodity:`SYN; contractYM:202006; expiry:2020.06.20; firstDate:2020.01.01; date:d; open:130#0f; high:130#0f; low:130#0f; settle:48f+0.04*til 130; volume:130#80);
futures:c1,c2;
.cfg.rolls[`SYN]:`type`rollDays`W`priceSource`method!(`days_before_expiry;10;5;`settle;`difference);

/ --- days_before_expiry: held -> window blend -> rolled ---
r19:.roll.active[2020.03.01;`SYN];          / 19 days to c1 expiry (> 10+5) -> held c1
chk[(r19`active)=202003; "held the nearest contract when far from its roll point"];
chk[(r19`reason)~`held; "reason is `held outside the window"];
r13:.roll.active[2020.03.07;`SYN];          / 13 days -> in (10,15] window: frac=(13-10)/5=0.6
chk[(r13`reason)~`rolling; "inside the roll window the reason is `rolling"];
chk[apx[0.6; (r13`weights) 202003; 1e-12]; "front weight = (heldD2e-rollDays)/W during the window"];
chk[apx[0.4; (r13`weights) 202006; 1e-12]; "next weight = 1 - front weight (the blend sums to 1)"];
r10:.roll.active[2020.03.10;`SYN];          / exactly 10 days = rollDays -> rolled out of c1
chk[(r10`active)=202006; "rolls exactly rollDays before expiry (out of the front at d2e<=rollDays)"];

/ --- THE KEY INVARIANT: a planted FUTURE volume spike does NOT change a past decision ---
.cfg.rolls[`VS]:`type`rollDays`W`priceSource`method!(`volume_switch;10;3;`settle;`difference);
vbase:update commodity:`VS from c1,c2;       / c1 vol 100 > c2 vol 80 -> front active
futures:vbase;
before:(.roll.active[2020.02.01;`VS])`active;
futures:vbase, update volume:999999 from select from vbase where contractYM=202006, date=2020.05.01;   / FUTURE spike
after:(.roll.active[2020.02.01;`VS])`active;
chk[before~after; "a FUTURE volume spike (date>asOf) must NOT change a past roll decision (as-of-only)"];
chk[before=202003; "as-of volume_switch holds the higher-as-of-volume front contract"];

/ --- volume_switch rolls when the next contract's as-of volume crosses the front's ---
c2hi:update volume:300 from c2;              / c2 now more liquid than c1 (100)
futures:update commodity:`VS from c1,c2hi;
chk[202006=(.roll.active[2020.02.01;`VS])`active; "volume_switch rolls to the higher-as-of-volume next contract"];

/ --- oi_switch is recognised but flagged unavailable (no OI column) ---
.cfg.rolls[`OI]:`type`rollDays`W`priceSource`method!(`oi_switch;10;5;`settle;`difference);
futures:update commodity:`OI from c1,c2;
chk[caught[{.roll.active[2020.02.01;`OI]};`]; "oi_switch must be flagged UNAVAILABLE (no open-interest data)"];

-1 "test_roll_active: PASS";
