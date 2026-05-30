\l core/init.q
/ loadContract on synthetic in-memory CSV text using the REAL header
/ (Time,Open,High,Low,Latest,Volume). No file dependency.
csvText:"Time,Open,High,Low,Latest,Volume",
    "\n2018-12-27T06:00:00+0000,48.3,48.3,47.34,47.52,4505",
    "\n2018-12-28T06:00:00+0000,47.94,47.94,47.81,47.82,1315",
    "\n2019-12-17T06:00:00+0000,60.1,60.5,59.9,60.4,1000",
    "\n2019-12-18T06:00:00+0000,60.55,61.18,60.32,60.93,131071",
    "\n2019-12-19T06:00:00+0000,60.86,61.47,60.79,61.22,20830";
tbl:.parser.crude.loadContract csvText;

.testutil.assertTrue[5=count tbl;"5 data rows (header dropped)"];
.testutil.assertTrue[(`date`open`high`low`settle`volume)~cols tbl;"stable OHLCV schema"];
/ Date parsing = first 10 chars of the ISO Time (offset/time ignored).
.testutil.assertTrue[2018.12.27=first tbl`date;"date from first 10 chars of Time"];
/ settle is the Latest column.
.testutil.assertNear[(first tbl`settle);47.52;1e-9;"settle = Latest"];
.testutil.assertNear[(last tbl`settle);61.22;1e-9;"last settle = Latest"];
/ expiry / firstDate derived from data (max / min date).
.testutil.assertTrue[2019.12.19=.parser.crude.contractExpiry tbl;"expiry = max date"];
.testutil.assertTrue[2018.12.27=.parser.crude.contractFirstDate tbl;"firstDate = min date"];
/ column types: settle float, volume long.
.testutil.assertTrue[9h=type tbl`settle;"settle is a float vector"];
.testutil.assertTrue[7h=type tbl`volume;"volume is a long vector"];
.testutil.assertTrue[20830=last tbl`volume;"volume parsed (last row)"];

/ Header-only CSV (with a newline) -> controlled error.
emptyErr:@[.parser.crude.loadContract;"Time,Open,High,Low,Latest,Volume\n";{`ERROR}];
.testutil.assertTrue[emptyErr~`ERROR;"header-only CSV rejected"];

-1 "PASS test_parser_crude_loadcontract: rows=",(string count tbl),", expiry=",string .parser.crude.contractExpiry tbl;
