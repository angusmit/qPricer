\l lib/init.q

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
rollResult:.backtest.commodity.longFrontRoll[futuresCurve;2024.01.02;2024.01.04;20;1f];
.testutil.assertTrue[(count rollResult)>0;"roll has rows"];
.testutil.assertTrue[`dailyPnl in key rollResult 0;"has dailyPnl"];
.testutil.assertTrue[`cumulativePnl in key rollResult 0;"has cumulativePnl"];
.testutil.assertTrue[`rolled in key rollResult 0;"has rolled flag"];

/ Check PnL: day2 pnl = 1000 * (73.0 - 72.5) = 500
secondRow:rollResult 1;
.testutil.assertNear[secondRow`dailyPnl;500f;1f;"day2 pnl ~500"];

-1 "PASS test_front_roll_backtest: rows=",string count rollResult;
