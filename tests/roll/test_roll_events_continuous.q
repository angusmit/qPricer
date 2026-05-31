\l core/init.q
/ ============================================================================
/ test_roll_events_continuous.q - roll events + the analytics-only continuous view. R11.
/ Synthetic. Roll events are emitted exactly at the active-contract changes (right date,
/ from/to); the continuous-adjusted series has NO gap at the roll (it uses the new contract's
/ own move, not the contract-switch jump); the rollEvents table is a single flat (unpartitioned)
/ table.
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};
apx:{[a;b;tol] tol>abs a-b};

d:2020.01.01+til 130;
c1:([] commodity:`SYN; contractYM:202003; expiry:2020.03.20; firstDate:2020.01.01; date:d; open:130#0f; high:130#0f; low:130#0f; settle:50f+0.05*til 130; volume:130#100);
c2:([] commodity:`SYN; contractYM:202006; expiry:2020.06.20; firstDate:2020.01.01; date:d; open:130#0f; high:130#0f; low:130#0f; settle:48f+0.04*til 130; volume:130#80);
futures:c1,c2;
.cfg.rolls[`SYN]:`type`rollDays`W`priceSource`method!(`days_before_expiry;10;5;`settle;`difference);

/ --- roll events: exactly one c1 -> c2 change, at the right date ---
ev:.roll.events[`SYN;2020.01.01;2020.06.30];
chk[1=count ev; "exactly one roll event for a 2-contract schedule (c1 -> c2)"];
chk[(first ev`fromContract)=202003; "the event rolls FROM the front contract"];
chk[(first ev`toContract)=202006; "the event rolls TO the next contract"];
/ the change date is where .roll.active first becomes c2.
firstC2:first (.data.hdb.dates `SYN) where 202006=({[d](.roll.active[d;`SYN])`active} each .data.hdb.dates `SYN);
chk[(first ev`date)=firstC2; "the event date is exactly where the active contract changes"];

/ --- the rollEvents table is a single flat (unpartitioned) table; persistable with p# ---
chk[98h=type ev; "rollEvents is a single flat (unpartitioned) table"];
chk[all `date`commodity`fromContract`toContract`reason in cols ev; "rollEvents has the expected columns"];

/ --- continuous-adjusted view: NO gap at the roll (uses the new contract's own move) ---
cont:.roll.continuous[`SYN];
ri:first where cont`rolled;
adjMove:(cont[`adjusted] ri)-cont[`adjusted] ri-1;
rawMove:(cont[`rawFront] ri)-cont[`rawFront] ri-1;
chk[apx[0.04; adjMove; 1e-9]; "the adjusted series moves by the NEW contract's own daily move (0.04), not the contract gap"];
chk[1f < abs rawMove; "the raw front series JUMPS at the roll (the ~2 contract gap) - which the adjustment removes"];
chk[(abs adjMove) < abs rawMove; "the adjusted move is far smaller than the raw contract-switch jump (no gap)"];
chk[(count cont)=count .data.hdb.dates `SYN; "the continuous series spans every date"];

-1 "test_roll_events_continuous: PASS";
