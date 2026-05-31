\l core/init.q
/ ============================================================================
/ test_template_registry.q - templates are first-class plug-ins via R2's `.registry.*`
/ (`template` kind). Round-trip register/get/list; a DEFICIENT template manifest is
/ REJECTED by the conformance check (the contract bites for templates too). Research OS R6.
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};

/ --- the two shipped templates are registered + conform ---
chk[`directionalSignal in .template.list[]; "directionalSignal must be registered"];
chk[`relativeValue in .template.list[]; "relativeValue must be registered"];
chk[.template.conforms `directionalSignal; "directionalSignal must conform to the template contract"];
chk[.template.conforms `relativeValue; "relativeValue must conform to the template contract"];

/ --- round-trip a synthetic template ---
fn:{[inputs] `pnl`validation`meta!(([] date:enlist 2020.01.01; pnl:enlist 0f); ()!(); ()!())};
.template.register[`synTmpl;fn;`contractVersion`shape`in`out!(1;`synthetic;()!();`pnl`validation`meta!`table`dict`dict)];
chk[`synTmpl in .template.list[]; "the synthetic template must appear in the list"];
chk[((.template.get[`synTmpl])`fn)~fn; "get must return the template fn handle"];
chk[.template.conforms `synTmpl; "a well-formed template manifest conforms"];

/ --- run dispatches to the registered fn ---
r:.template.run[`synTmpl; ()!()];
chk[(`pnl`validation`meta)~key r; "run must return pnl/validation/meta"];

/ --- a DEFICIENT manifest (missing the required output key `pnl) is REJECTED ---
.template.register[`badTmpl;fn;`contractVersion`shape`in`out!(1;`synthetic;()!();`validation`meta!`dict`dict)];
chk[not .template.conforms `badTmpl; "a template manifest missing `pnl in out must be rejected"];
chk[not (first exec pass from .contracts.verify[] where kind=`template, name=`badTmpl); ".contracts.verify must FLAG the deficient template"];

/ cleanup: drop the synthetic templates so .contracts.verify[] stays clean for other tests.
.registry.unregister[`template;`synTmpl];
.registry.unregister[`template;`badTmpl];
chk[all (.contracts.verify[])`pass; "after cleanup, every registered plug-in still conforms"];

-1 "test_template_registry: PASS";
