\l core/init.q
/ curveAt on a synthetic long table (3 contracts, differing expiries). Built
/ column-wise in q; no CSV files. asof = 2020.01.06, all three alive.
asofDate:2020.01.06;
/ contract 202001 carries a SECOND row on 2020.01.07 (settle 99) to prove
/ curveAt picks the asof settle (50), not some other date's settle.
contractYM:202001 202002 202003 202001;
expiry:    2020.01.20 2020.02.20 2020.03.20 2020.01.20;
firstDate: 4#2019.12.01;
date:      2020.01.06 2020.01.06 2020.01.06 2020.01.07;
settle:    50.0 51.0 52.0 99.0;
longTable:flip `contractYM`expiry`firstDate`date`settle`open`high`low`volume!(
    contractYM;expiry;firstDate;date;settle;settle;settle;settle;4#1000);

curve:.parser.crude.curveAt[longTable;asofDate];
.testutil.assertTrue[3=count curve;"3 alive+quoted contracts on asof"];
.testutil.assertTrue[(`tenor`price`contractYM`expiry)~cols curve;"curve (tenor,price,...) schema"];
.testutil.assertTrue[(curve`tenor)~asc curve`tenor;"tenors ascending"];
.testutil.assertTrue[(202001 202002 202003)~curve`contractYM;"contracts ordered by ascending tenor"];

/ tenor arithmetic: 202001 expiry 2020.01.20, asof 2020.01.06 -> 14/365 years.
.testutil.assertNear[(first curve`tenor);14%365f;1e-9;"tenor = (expiry-asof)/365"];
/ 202003 expiry 2020.03.20 -> (2020.03.20 - 2020.01.06) = 74 days.
.testutil.assertNear[(last curve`tenor);74%365f;1e-9;"far tenor = 74/365"];
/ price = settle ON asof (50), not the 2020.01.07 row's 99.
.testutil.assertNear[(first curve`price);50.0;1e-9;"price = asof settle, not later date"];
.testutil.assertTrue[(50.0 51.0 52.0)~curve`price;"prices in tenor order"];

/ A contract alive but with NO quote on asof is dropped (not crashing): drop the
/ 202002 asof row, leaving it quote-less on asof.
sparseTable:longTable where not (longTable[`contractYM]=202002) and longTable[`date]=asofDate;
sparseCurve:.parser.crude.curveAt[sparseTable;asofDate];
.testutil.assertTrue[2=count sparseCurve;"contract with no asof quote dropped"];
.testutil.assertTrue[not 202002 in sparseCurve`contractYM;"202002 absent (no asof settle)"];

-1 "PASS test_parser_crude_curve: tenors=",(-3!curve`tenor),", prices=",-3!curve`price;
