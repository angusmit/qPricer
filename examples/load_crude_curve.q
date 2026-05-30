\l lib/init.q
/ ============================================================================
/ load_crude_curve.q - real-data example for the v0.50 CRUDE futures parser.
/ ----------------------------------------------------------------------------
/ Reads the REAL Barchart WTI daily settle CSVs under data/barchart/CRUDE/ and
/ prints a forward-curve snapshot (tenor, price) for a sample as-of date, with a
/ contango/backwardation note.
/ ----------------------------------------------------------------------------
/ This is an EXAMPLE, not a test. The CSVs are USER-SUPPLIED and gitignored
/ (see .gitignore: /data/barchart/CRUDE/*.csv) - they are NOT committed. The
/ test suite exercises the parser on synthetic in-memory CSV text instead, so it
/ does not depend on this folder.
/ ============================================================================

crudeDir:"data/barchart/CRUDE";
-1 "Loading CRUDE futures curve from ",crudeDir," ...";

longTable:.parser.crude.loadAll crudeDir;
-1 "Loaded ",(string count longTable)," daily settle rows across ",
    (string count distinct longTable`contractYM)," contracts.";
-1 "Date range: ",(string min longTable`date)," .. ",string max longTable`date;

/ Pick a sample as-of date that sits inside the data range and has several alive
/ contracts (early 2020). Fall back to the median traded date if absent.
sampleAsof:2020.01.06;
tradedDates:asc distinct longTable`date;
if[not sampleAsof in tradedDates;
   sampleAsof:tradedDates (count tradedDates) div 2];

-1 "";
-1 "Forward curve as of ",(string sampleAsof),":";
curve:.parser.crude.curveAt[longTable;sampleAsof];
show curve;

/ Contango (upward) vs backwardation (downward): compare the near vs far price.
nearPrice:first curve`price;
farPrice:last curve`price;
shape:$[farPrice>nearPrice;"CONTANGO (far > near)";
        farPrice<nearPrice;"BACKWARDATION (far < near)";
        "FLAT"];
-1 "";
-1 "Near tenor ",(string first curve`tenor)," yr  price ",string nearPrice;
-1 "Far  tenor ",(string last curve`tenor)," yr  price ",string farPrice;
-1 "Curve shape: ",shape;
