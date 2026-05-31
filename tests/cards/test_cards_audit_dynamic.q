\l core/init.q
/ ============================================================================
/ test_cards_audit_dynamic.q - the complete audit (.cards.audit / .cards.__registered).
/ v0.71 wiring: the audit walks ALL kinds in the R2 registry dynamically, so a NEW kind
/ (the R6 `template` kind, or any future kind) can never be silently missed. Synthetic.
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};

/ --- the live audit already covers the R6 `template` kind (not a hardcoded list) ---
a0:.cards.audit[];
chk[`template in a0`kind; "the audit must dynamically cover the R6 template kind"];
chk[all `model`signal`calibrator`fillModel`template`strategy in a0`kind; "the audit covers every registry kind + strategy"];

/ --- register a brand-new synthetic kind: the audit must pick it up with NO code change ---
.registry.new[`synKind; `version`requiredIn`requiredOut!(1; ()!(); (enlist `x)!enlist `float)];
.registry.register[`synKind;`synCap;{x};`contractVersion`in`out!(1; ()!(); (enlist `x)!enlist `float)];
a:.cards.audit[];
chk[`synKind in a`kind; "a newly-registered kind must appear in the audit (dynamic walk)"];
chk[`synCap in exec name from a where kind=`synKind; "the new kind's capability is covered"];
chk[not first exec pass from a where kind=`synKind, name=`synCap; "an uncarded new capability is flagged (the audit bites)"];

/ cleanup: drop the synthetic kind so the registry stays clean for other tests.
.registry.dropKind[`synKind];
chk[not `synKind in (.cards.audit[])`kind; "after dropKind the synthetic kind is gone from the audit"];

-1 "test_cards_audit_dynamic: PASS";
