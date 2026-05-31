\l core/init.q
/ ============================================================================
/ test_contracts_conform.q - the conformance check that CAN FAIL (.registry.conforms
/ / .contracts.verify). Synthetic kind, no HDB. Research OS R2.
/ A CONFORMING manifest passes; a DEFICIENT manifest (missing a required output key)
/ is REJECTED; a TYPE-TAG mismatch is rejected; a wrong contract version is rejected;
/ .contracts.verify[] flags the bad one and passes the good one. (A check that always
/ passes is theatre - this proves it bites.)
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};

/ a synthetic kind: requires in a,b (floats) and out x (float), y (symbol).
.registry.new[`tconf; `version`requiredIn`requiredOut!(1; `a`b!`float`float; `x`y!`float`symbol)];

/ --- a conforming manifest (declares all required keys/types; extra `c is fine) ---
.registry.register[`tconf;`good;{x};
    `contractVersion`in`out!(1; `a`b`c!`float`float`float; `x`y!`float`symbol)];
chk[.registry.conforms[`tconf;`good]; "a conforming manifest must pass"];

/ --- DEFICIENT: missing required output key `y -> REJECTED (the contract bites) ---
.registry.register[`tconf;`missingOut;{x};
    `contractVersion`in`out!(1; `a`b!`float`float; (enlist `x)!enlist `float)];
chk[not .registry.conforms[`tconf;`missingOut]; "a manifest missing a required OUTPUT key must be rejected"];

/ --- DEFICIENT: missing required input key `b -> REJECTED ---
.registry.register[`tconf;`missingIn;{x};
    `contractVersion`in`out!(1; (enlist `a)!enlist `float; `x`y!`float`symbol)];
chk[not .registry.conforms[`tconf;`missingIn]; "a manifest missing a required INPUT key must be rejected"];

/ --- TYPE MISMATCH: out `y declared float but contract requires symbol -> REJECTED ---
.registry.register[`tconf;`badType;{x};
    `contractVersion`in`out!(1; `a`b!`float`float; `x`y!`float`float)];
chk[not .registry.conforms[`tconf;`badType]; "a type-tag mismatch on a required key must be rejected"];

/ --- WRONG CONTRACT VERSION -> REJECTED ---
.registry.register[`tconf;`badVer;{x};
    `contractVersion`in`out!(2; `a`b!`float`float; `x`y!`float`symbol)];
chk[not .registry.conforms[`tconf;`badVer]; "an incompatible contract version must be rejected"];

/ --- .contracts.verify[] flags exactly the bad ones and passes the good one ---
v:.contracts.verify[];
vt:select from v where kind=`tconf;
chk[(first exec pass from vt where name=`good); "verify must pass the conforming plug-in"];
chk[not any exec pass from vt where name in `missingOut`missingIn`badType`badVer; "verify must FAIL every deficient plug-in"];
chk[4=count select from vt where not pass; "verify must flag all four deficient plug-ins"];

/ cleanup: drop the synthetic kind so .contracts.verify[] stays clean for other tests.
.registry.dropKind[`tconf];
chk[all (.contracts.verify[])`pass; "after cleanup, every remaining (real) capability still conforms"];

-1 "test_contracts_conform: PASS";
