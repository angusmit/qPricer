\l lib/init.q
/ Generalised futures parser: parses any product tag from the Day_<TAG>_YYYYMM00
/ filename, the tag-independent CSV path is shared with crude (byte-identical),
/ and the crude parser is unchanged for its own tag. Synthetic CSV text only.
.testutil.assertTrue["Day_GAS_"~.parser.futures.__fnPrefix[`GAS];"futures filename prefix for a tag"];
.testutil.assertTrue["Day_HO_"~.parser.futures.__fnPrefix[`HO];"futures filename prefix for another tag"];

/ Filename month parsing works for a non-crude tag.
.testutil.assertTrue[(2020 1)~.parser.futures.contractMonthFromFilename "Day_GAS_20200100.csv";"GAS filename -> (2020;1)"];
.testutil.assertTrue[(2021 12)~.parser.futures.contractMonthFromFilename "data/barchart/GAS/Day_GAS_20211200.csv";"GAS path -> (2021;12)"];
/ Crude filename parsing is unchanged.
.testutil.assertTrue[(2020 1)~.parser.crude.contractMonthFromFilename "Day_CRUDE_20200100.csv";"crude filename parsing unchanged"];

/ The CSV-content parse is tag-independent: futures.loadContract == crude.loadContract.
csvText:"Time,Open,High,Low,Latest,Volume",
    "\n2020-01-02T06:00:00+0000,2.10,2.20,2.00,2.15,1000",
    "\n2020-01-03T06:00:00+0000,2.15,2.25,2.10,2.20,1200",
    "\n2020-01-06T06:00:00+0000,2.20,2.30,2.15,2.25,1500";
futuresTbl:.parser.futures.loadContract csvText;
crudeTbl:.parser.crude.loadContract csvText;
.testutil.assertTrue[futuresTbl~crudeTbl;"futures.loadContract byte-identical to crude.loadContract"];
.testutil.assertTrue[3=count futuresTbl;"3 parsed rows"];
.testutil.assertNear[(first futuresTbl`settle);2.15;1e-9;"settle from Latest column"];

-1 "PASS test_parser_futures";
