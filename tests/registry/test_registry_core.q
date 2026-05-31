\l core/init.q
/ ============================================================================
/ test_registry_core.q - the generic registry primitive (.registry.*). Synthetic
/ kind, no HDB. Research OS R2. Round-trip register/get/list + defined duplicate-name
/ (overwrite-with-warn, latest wins) and unknown-name / unknown-kind (signal) behaviour.
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};
caught:{[f;arg] @[f;arg;{[e] 1b}]~1b};   / true iff f arg signalled

/ a synthetic kind with a trivial contract (dropped at the end for test independence).
.registry.new[`tcore; `version`requiredIn`requiredOut!(1; (enlist `a)!enlist `float; (enlist `x)!enlist `float)];

/ --- round-trip: register -> get -> list ---
fnA:{x+1};
.registry.register[`tcore;`plugA;fnA;`contractVersion`in`out!(1;(enlist `a)!enlist `float;(enlist `x)!enlist `float)];
chk[(.registry.get[`tcore;`plugA]`fn)~fnA; "get must return the registered fn handle"];
chk[`plugA in .registry.list `tcore; "list must include the registered name"];
chk[1=count .registry.list `tcore; "exactly one plug-in registered so far"];

/ --- duplicate name: overwrite-with-warn, the LATEST registration wins (no dup row) ---
fnA2:{x+2};
.registry.register[`tcore;`plugA;fnA2;`contractVersion`in`out!(1;(enlist `a)!enlist `float;(enlist `x)!enlist `float)];
chk[1=count .registry.list `tcore; "a duplicate name must NOT add a second entry"];
chk[(.registry.get[`tcore;`plugA]`fn)~fnA2; "the latest registration must win"];

/ --- unknown name -> signal ---
chk[caught[.registry.get[`tcore;];`nope]; "get of an unknown name must signal"];

/ --- unknown kind -> signal (get / list / register) ---
chk[caught[.registry.list;`nokind]; "list of an unknown kind must signal"];
chk[caught[.registry.get[`nokind;];`x]; "get on an unknown kind must signal"];

/ --- a second plug-in, then unregister ---
.registry.register[`tcore;`plugB;{x};`contractVersion`in`out!(1;(enlist `a)!enlist `float;(enlist `x)!enlist `float)];
chk[2=count .registry.list `tcore; "two plug-ins registered"];
.registry.unregister[`tcore;`plugA];
chk[(enlist `plugB)~.registry.list `tcore; "unregister must remove exactly the named plug-in"];

/ cleanup: drop the synthetic kind so the global registry stays clean for other tests.
.registry.dropKind[`tcore];
chk[not `tcore in .registry.kinds[]; "dropKind must remove the synthetic kind"];

-1 "test_registry_core: PASS";
