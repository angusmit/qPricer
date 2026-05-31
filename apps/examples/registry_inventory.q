\l core/init.q
/ ============================================================================
/ registry_inventory.q - the plug-in spine, made visible (Research OS R2). Lists each
/ kind's registered capabilities + runs .contracts.verify[], printing the inventory +
/ the conformance result: a one-screen view that the spine is populated and every
/ registered capability conforms to its declared contract. No HDB needed. exit 0;.
/ ============================================================================
-1 "qFDM capability registry - the plug-in spine";
-1 "================================================";
-1 "";
{[kind]
    c:.registry.contract kind;
    -1 (upper string kind)," (contract v",(string c`version),"):";
    -1 "  required in : ",(", " sv string key c`requiredIn);
    -1 "  required out: ",(", " sv string key c`requiredOut);
    names:.registry.list kind;
    {[kind;nm] m:(.registry.get[kind;nm])`manifest;
        -1 "    - ",(string nm),"  ",$[`description in key m; m`description; ""]}[kind] each names;
    -1 ""} each .registry.kinds[];

-1 "Conformance (.contracts.verify[]):";
v:.contracts.verify[];
show select kind,name,pass,reason from v;
-1 "";
-1 $[all v`pass;
    "ALL ",(string count v)," registered capabilities conform to their contracts.";
    "WARNING: ",(string count select from v where not pass)," capability(ies) do NOT conform - see above."];
-1 "";
-1 "(The spine is one generic registry + a versioned contract per kind + a conformance";
-1 " check that can REFUSE a non-conforming plug-in. R5/R6/R7 snap in by declaring a contract.)";
exit 0;
