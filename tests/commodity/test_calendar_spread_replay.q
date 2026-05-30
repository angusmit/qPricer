\l core/init.q

/ Synthetic WTI futures: 3 dates, 3 contracts per date
futuresTable:();
dates:2024.01.02 2024.01.03 2024.01.04;
contracts:`CLF24`CLG24`CLH24;
deliveryMonths:2024.01.15 2024.02.15 2024.03.15;
expiryDates:2024.01.18 2024.02.16 2024.03.15;
basePrices:(72.5 73.0 73.8;73.0 73.5 74.2;71.5 72.0 72.7);
dIdx:0;
while[dIdx<3;
    cIdx:0;
    while[cIdx<3;
        futuresTable:futuresTable,enlist `snapshotDate`commodity`exchange`contractCode`deliveryMonth`expiryDate`settlementPrice`volume`openInterest`contractMultiplier!(
            dates dIdx;`WTI;`NYMEX;contracts cIdx;deliveryMonths cIdx;expiryDates cIdx;basePrices[dIdx;cIdx];50000f;200000f;1000f);
        cIdx+:1];
    dIdx+:1];

futuresCurve:.commodity.futures.buildCurveByDate futuresTable;
spreadResult:.backtest.commodity.calendarSpreadReplay[futuresCurve;1;3;2024.01.02;2024.01.04;1f];
.testutil.assertTrue[(count spreadResult)>0;"spread has rows"];
firstRow:spreadResult 0;
.testutil.assertTrue[`spread in key firstRow;"has spread"];
.testutil.assertNear[firstRow`spread;72.5-73.8;0.01;"spread = near - far"];

/ Cumulative works
lastRow:spreadResult (count spreadResult)-1;
.testutil.assertTrue[not null lastRow`cumulativePnl;"cumPnl populated"];

-1 "PASS test_calendar_spread_replay: rows=",string count spreadResult;
