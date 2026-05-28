\l lib/init.q
/ fromBarchart maps an already-normalised in-memory Barchart contract series to the
/ standard path schema. No file IO; we synthesise a minimal normalised table.

synthNormalised:([] snapshotDate:2024.01.02 2024.01.03 2024.01.04 2024.01.05;
                    underlying:4#`aapl;
                    contractId:4#`aapl_C100;
                    expiryDate:4#2024.02.16;
                    optionType:4#`call;
                    strike:4#100f;
                    spot:185.0 186.5 184.2 187.1;
                    marketPrice:5.5 6.1 5.4 6.2;
                    impliedVolatility:0.28 0.27 0.29 0.275;
                    status:4#`OK);

extraContract:([] snapshotDate:enlist 2024.01.02;
                   underlying:enlist `aapl;
                   contractId:enlist `aapl_P95;
                   expiryDate:enlist 2024.02.16;
                   optionType:enlist `put;
                   strike:enlist 95f;
                   spot:enlist 185.0;
                   marketPrice:enlist 1.2;
                   impliedVolatility:enlist 0.31;
                   status:enlist `OK);

mixed:synthNormalised,extraContract;

pathTbl:.strategy.path.fromBarchart[mixed;`aapl_C100];

requiredCols:`stepIndex`stepDate`spot`volatility`riskFreeRate`dividendYield`marketPrice`status;
.testutil.assertTableColumns[pathTbl;requiredCols;"barchart path schema"];
.testutil.assertTrue[4=count pathTbl;"4 rows from 4 contract snapshots"];
.testutil.assertTrue[(til 4)~pathTbl`stepIndex;"stepIndex = 0..3"];

datesAsc:asc synthNormalised`snapshotDate;
.testutil.assertTrue[(datesAsc)~pathTbl`stepDate;"stepDate sorted ascending by snapshotDate"];
.testutil.assertTrue[185.0=first pathTbl`spot;"first spot = first snapshot spot"];
.testutil.assertTrue[(pathTbl`volatility)~synthNormalised`impliedVolatility;"volatility = impliedVolatility"];
.testutil.assertTrue[(pathTbl`marketPrice)~synthNormalised`marketPrice;"marketPrice carried through"];

unknownResult:.[.strategy.path.fromBarchart;(mixed;`aapl_XYZ);{`ERROR}];
.testutil.assertTrue[unknownResult~`ERROR;"unknown contractId rejected"];

.strategy.path.validate pathTbl;
-1 "  validate passed on barchart-adapted path";

-1 "PASS test_strategy_path_barchart_adapter: rows=",string[count pathTbl];
