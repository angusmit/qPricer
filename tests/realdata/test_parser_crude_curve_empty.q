\l core/init.q
/ curveAt controlled-error cases: asof before any contract is alive and asof
/ after all contracts have expired. Synthetic long table, no files.
contractYM:202001 202002;
expiry:    2020.01.20 2020.02.20;
firstDate: 2#2019.12.01;
date:      2020.01.06 2020.01.06;
settle:    50.0 51.0;
longTable:flip `contractYM`expiry`firstDate`date`settle!(
    contractYM;expiry;firstDate;date;settle);

/ Before any firstDate -> nothing alive -> controlled error.
beforeErr:@[.parser.crude.curveAt[longTable;];2019.10.01;{`ERROR}];
.testutil.assertTrue[beforeErr~`ERROR;"asof before all contracts -> controlled error"];

/ After all expiries -> nothing alive -> controlled error.
afterErr:@[.parser.crude.curveAt[longTable;];2020.06.01;{`ERROR}];
.testutil.assertTrue[afterErr~`ERROR;"asof after all expiries -> controlled error"];

/ curveHistory over a mix of valid + invalid asof dates keeps only the valid
/ curve (invalid dates skipped, not fatal).
hist:.parser.crude.curveHistory[longTable;2019.10.01 2020.01.06 2020.06.01];
.testutil.assertTrue[(`asofDate`tenor`price`contractYM`expiry)~cols hist;"curveHistory schema (asofDate, not asof)"];
.testutil.assertTrue[all 2020.01.06=hist`asofDate;"only the one valid asof survives"];
.testutil.assertTrue[2=count hist;"2 alive contracts on the single valid asof"];

/ All-invalid asof dates -> controlled error.
allBadErr:@[.parser.crude.curveHistory[longTable;];2019.10.01 2020.06.01;{`ERROR}];
.testutil.assertTrue[allBadErr~`ERROR;"curveHistory with no valid dates -> error"];

-1 "PASS test_parser_crude_curve_empty";
