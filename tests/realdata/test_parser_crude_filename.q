\l core/init.q
/ Parse delivery year/month from the Barchart CRUDE filename convention
/ Day_CRUDE_YYYYMM00.csv (trailing 00 is Barchart's day filler).
ym:.parser.crude.contractMonthFromFilename "Day_CRUDE_20200100.csv";
.testutil.assertTrue[ym~2020 1;"Day_CRUDE_20200100.csv -> (2020;1)"];
ym2:.parser.crude.contractMonthFromFilename "data/barchart/CRUDE/Day_CRUDE_20211200.csv";
.testutil.assertTrue[ym2~2021 12;"path with dir + Dec -> (2021;12)"];

/ CL month code mapping (F G H J K M N Q U V X Z).
.testutil.assertTrue["F"=.parser.crude.monthCode 1;"month 1 -> F"];
.testutil.assertTrue["Z"=.parser.crude.monthCode 12;"month 12 -> Z"];

/ Invalid month in filename -> controlled error.
badMonth:@[.parser.crude.contractMonthFromFilename;"Day_CRUDE_20201300.csv";{`ERROR}];
.testutil.assertTrue[badMonth~`ERROR;"month 13 rejected"];
/ Malformed filename -> controlled error.
badName:@[.parser.crude.contractMonthFromFilename;"random.csv";{`ERROR}];
.testutil.assertTrue[badName~`ERROR;"malformed filename rejected"];

-1 "PASS test_parser_crude_filename: ",(-3!ym),", ",-3!ym2;
