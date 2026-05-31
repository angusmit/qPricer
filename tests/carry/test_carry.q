\l core/init.q
/ ============================================================================
/ test_carry.q - carry / storage economics (.carry.*). R14. SYNTHETIC.
/ Implied carry is positive on a backwardated curve / negative on contango; convenience yield =
/ (r+storage) - carry on known inputs; cash-and-carry fair value is the known-answer; the carry signal
/ flips sign with the curve; the inventory-tightness proxy tracks the backwardation degree.
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};
apx:{[a;b;tol] tol>abs a-b};

/ a 6-contract curve snapshot as of 2025.01.05 with the given prices (front..CL6), descending=backwardation.
mexp:{"D"$(string x div 100),".",(-2#"0",string x mod 100),".20"};
mkCurve:{[prices] yms:202502+til 6;                       / 6 monthly contracts 202502..202507
    ds:2025.01.01+til 5;                                  / a few dates up to the asOf (2025.01.05)
    raze {[ym;px;ds;mexp] e:mexp ym; n:count ds;
        flip `commodity`contractYM`expiry`firstDate`date`settle`open`high`low`volume!(
            (n#`SYN);(n#ym);(n#e);(n#first ds);ds;(n#px);(n#px);(n#px);(n#px);(n#5000))}[;;ds;mexp]'[yms;prices]};
asOf:2025.01.05;
r:.cfg.carry`riskFreeRate; storage:.cfg.carry`storageCost; coc:r+storage;

/ --- backwardation: prices descending (front highest) ---
futures:mkCurve 80 79 78 77 76 75f;
b:.carry.features[asOf;`SYN];
chk[`backwardation=b`classification; "descending curve classified backwardation"];
chk[0f<b`impliedCarry; "implied carry is POSITIVE on a backwardated curve"];
chk[0f<b`inventoryTightness; "inventory-tightness proxy is positive (tight) in backwardation"];
/ formula known-answers (using the output's own front/deferred/tenorGap)
chk[apx[(log (b`front)%b`deferred)%b`tenorGap; b`impliedCarry; 1e-9]; "implied carry == ln(front/deferred)/T"];
chk[apx[coc-b`impliedCarry; b`convenienceYield; 1e-12]; "convenience yield == (r+storage) - implied carry"];
chk[apx[(b`front)*exp coc*b`tenorGap; b`fairValueCarry; 1e-9]; "cash-and-carry fair value == front*exp((r+storage)*T)"];
chk[apx[(b`impliedCarry)-coc; b`carrySignal; 1e-12]; "carry signal == implied carry - cost of carry"];
chk[(`convenienceYield`fairValueCarry)~b`assumptionDependent; "convenience yield + fair value flagged assumption-dependent"];

/ --- contango: prices ascending (front lowest) ---
futures:mkCurve 70 71 72 73 74 75f;
c:.carry.features[asOf;`SYN];
chk[`contango=c`classification; "ascending curve classified contango"];
chk[0f>c`impliedCarry; "implied carry is NEGATIVE on a contango curve"];
chk[0f>c`carrySignal; "carry signal is negative in contango (flips sign with the curve)"];

/ --- the signals flip / track the curve across the two regimes ---
chk[(b`carrySignal)>c`carrySignal; "the carry signal flips: higher in backwardation than contango"];
chk[(b`inventoryTightness)>c`inventoryTightness; "inventory-tightness is higher (tighter) in backwardation"];

-1 "test_carry: PASS";
