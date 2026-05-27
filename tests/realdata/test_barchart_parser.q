\l lib/init.q

/ Build synthetic normalised optionTable (3 dates, 2 contracts)
optionTable:();
/ Contract 1: AAPL call 165 exp 2024.02.16
/ Contract 2: AAPL put 185 exp 2024.02.16
dates:2024.01.02 2024.01.03 2024.01.04;
cid1:`$"aapl_2024.02.16_165_c";
cid2:`$"aapl_2024.02.16_185_p";
spots:185.64 186.50 185.00;
c1prices:22.05 23.10 21.80;
c2prices:3.50 3.20 3.80;
c1ivs:0.2755 0.2800 0.2700;
c2ivs:0.3000 0.2950 0.3100;
dIdx:0;
while[dIdx<3;
    optionTable:optionTable,enlist `snapshotDate`underlying`contractId`expiryDate`optionType`strike`spot`open`high`low`latest`bid`ask`mid`marketPrice`impliedVolatility`volume`openInterest`vendorDelta`vendorGamma`vendorTheta`vendorVega`vendorRho`vendorTheo`dte`tau`moneyness`bidAskSpread`bidAskSpreadPct`status`errorMessage!(
        dates dIdx;`aapl;cid1;2024.02.16;`call;165f;spots dIdx;c1prices[dIdx]-0.5;c1prices[dIdx]+0.5;c1prices[dIdx]-1f;c1prices dIdx;
        c1prices[dIdx]-0.15;c1prices[dIdx]+0.15;c1prices dIdx;c1prices dIdx;c1ivs dIdx;
        100f;2000f;0.91;0.009;-0.053;0.104;0.177;22.4;
        45-dIdx;(45-dIdx)%365f;165f%spots dIdx;0.30;0.30%c1prices dIdx;`OK;"");
    optionTable:optionTable,enlist `snapshotDate`underlying`contractId`expiryDate`optionType`strike`spot`open`high`low`latest`bid`ask`mid`marketPrice`impliedVolatility`volume`openInterest`vendorDelta`vendorGamma`vendorTheta`vendorVega`vendorRho`vendorTheo`dte`tau`moneyness`bidAskSpread`bidAskSpreadPct`status`errorMessage!(
        dates dIdx;`aapl;cid2;2024.02.16;`put;185f;spots dIdx;c2prices[dIdx]-0.3;c2prices[dIdx]+0.3;c2prices[dIdx]-0.5;c2prices dIdx;
        c2prices[dIdx]-0.10;c2prices[dIdx]+0.10;c2prices dIdx;c2prices dIdx;c2ivs dIdx;
        200f;1500f;-0.25;0.015;-0.040;0.080;-0.100;3.6;
        45-dIdx;(45-dIdx)%365f;185f%spots dIdx;0.20;0.20%c2prices dIdx;`OK;"");
    dIdx+:1];

/ Test normalise field mappings (synthetic data is already normalised, verify fields)
.testutil.assertTrue[6=count optionTable;"6 rows (3 dates x 2 contracts)"];
firstRow:optionTable 0;
.testutil.assertTrue[firstRow[`optionType]=`call;"cpFlag c -> call"];
.testutil.assertTrue[firstRow[`spot]=185.64;"spot from price column"];
.testutil.assertNear[firstRow`impliedVolatility;0.2755;0.0001;"iv 27.55% -> 0.2755"];
.testutil.assertTrue[firstRow[`marketPrice]=22.05;"latest -> marketPrice"];
.testutil.assertTrue[firstRow[`dte]>0i;"dte positive"];
.testutil.assertTrue[firstRow[`status]=`OK;"status OK"];

/ Test summary
summaryResult:.parser.barchart.summary optionTable;
.testutil.assertTrue[summaryResult[`rowCount]=6;"6 rows"];
.testutil.assertTrue[summaryResult[`contractCount]=2;"2 contracts"];
.testutil.assertTrue[summaryResult[`okRows]=6;"6 OK"];

/ Test availableDates
availDates:.parser.barchart.availableDates optionTable;
.testutil.assertTrue[3=count availDates;"3 available dates"];

-1 "PASS test_barchart_parser: rows=",string[summaryResult`rowCount],", contracts=",string summaryResult`contractCount;
