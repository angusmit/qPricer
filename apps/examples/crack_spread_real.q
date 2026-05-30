\l core/init.q
/ ============================================================================
/ crack_spread_real.q - real-data example (NOT a test), DATA-CONDITIONAL.
/ ----------------------------------------------------------------------------
/ IF crude + refined-product (RB gasoline, HO heating oil) files are present,
/ prices a crack-spread option on the real front curves with the EXISTING Kirk
/ spread pricer, using cross-correlation estimated from the real return series.
/ IF ABSENT, skips with a note. SPARK SPREAD IS OUT OF SCOPE (no power data).
/ CSVs are user-supplied and gitignored.
/ ============================================================================
crudeDir:"data/barchart/CRUDE";
rbDir:"data/barchart/RB";
hoDir:"data/barchart/HO";
hasDir:{0<count @[{key hsym `$x};x;{[e] ()}]};
if[not (hasDir rbDir) and hasDir hoDir;
    -1 "Crack-spread real demo SKIPPED - need data/barchart/RB and data/barchart/HO";
    -1 "(refined-product futures), which are not present (data-conditional).";
    -1 "The crackSpread strategy + spread pricers (Margrabe/Kirk/MC) are covered by";
    -1 "synthetic tests; drop Day_RB_*/Day_HO_* files in to run the real crack.";
    exit 0];

/ --- real branch (runs only when products data is present) ---
crudeLong:.parser.futures.loadAll[crudeDir;`CRUDE];
productLong:.parser.futures.loadAll[rbDir;`RB];
commonDates:asc distinct (distinct crudeLong`date) inter distinct productLong`date;
crudeCurve:.parser.futures.curveAt[crudeLong;last commonDates];
productCurve:.parser.futures.curveAt[productLong;last commonDates];
crudeFront:first crudeCurve`price;
productFront:first productCurve`price;
/ Estimate vols + correlation from overlapping front log-returns.
crudeFrontSeries:`float$exec first price by date from crudeLong where date in commonDates;
productFrontSeries:`float$exec first price by date from productLong where date in commonDates;
crudeRet:1_ deltas log crudeFrontSeries;
productRet:1_ deltas log productFrontSeries;
volC:(dev crudeRet)*sqrt 252f;
volP:(dev productRet)*sqrt 252f;
rho:crudeRet cor productRet;
/ Price a 1:1 crack-spread call (product - crude) with the existing Kirk pricer.
crackParams:`fwd1`fwd2`strike`expiry`vol1`vol2`correlation`riskFreeRate!(
    productFront;crudeFront;0f;0.25;volP;volC;rho;0.02);
crackPx:.commodity.spread.kirkPrice[`call;crackParams];
-1 "Real crack snapshot as of ",string last commonDates;
-1 "  product front=",(string productFront)," crude front=",(string crudeFront)," realized crack=",string productFront-crudeFront;
-1 "  vols (product/crude)=",(string volP),"/",(string volC)," correlation=",string rho;
-1 "  Kirk crack-spread call (K=0, 3m)=",string crackPx;
exit 0;
