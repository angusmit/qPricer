\l core/init.q

-1 "=============================================================================";
-1 " qFDM v0.32 Commodity Model Foundation Demo";
-1 "=============================================================================";
-1 "";

/ Build synthetic WTI futures curve
futuresTable:();
dates:2024.01.02+til 10;
contracts:`CLF24`CLG24`CLH24;
deliveryMonths:2024.01.15 2024.02.15 2024.03.15;
expiryDates:2024.01.18 2024.02.16 2024.03.15;
dIdx:0;
while[dIdx<10;
    basePx:72f+0.3*dIdx;
    cIdx:0;
    while[cIdx<3;
        futuresTable:futuresTable,enlist `snapshotDate`commodity`exchange`contractCode`deliveryMonth`expiryDate`settlementPrice`volume`openInterest`contractMultiplier!(
            dates dIdx;`WTI;`NYMEX;contracts cIdx;deliveryMonths cIdx;expiryDates cIdx;basePx+(0.5*cIdx)-(0.1*dIdx);50000f;200000f;1000f);
        cIdx+:1];
    dIdx+:1];

futuresCurve:.commodity.futures.buildCurveByDate futuresTable;

-1 "--- Futures Curve (day 1) ---";
day1:futuresCurve where (futuresCurve`snapshotDate)=2024.01.02;
fcIdx:0;
while[fcIdx<count day1;
    rw:day1 fcIdx;
    -1 "  ",string[rw`contractCode]," T",string[rw`tenorRank]," $",string[rw`settlementPrice]," DTE=",string rw`daysToExpiry;
    fcIdx+:1];
-1 "";

/ Front roll backtest
rollResult:.backtest.commodity.longFrontRoll[futuresCurve;2024.01.02;2024.01.11;15;1f];
-1 "--- Front Roll Backtest ---";
rrIdx:0;
while[rrIdx<count rollResult;
    rw:rollResult rrIdx;
    -1 "  ",string[rw`snapshotDate]," ",string[rw`currentContract]," ",string[rw`action]," pnl=",string[rw`dailyPnl]," cum=",string rw`cumulativePnl;
    rrIdx+:1];
-1 "";

/ Calendar spread replay
spreadResult:.backtest.commodity.calendarSpreadReplay[futuresCurve;1;3;2024.01.02;2024.01.11;1f];
-1 "--- Calendar Spread Replay ---";
srIdx:0;
while[srIdx<count spreadResult;
    rw:spreadResult srIdx;
    -1 "  ",string[rw`snapshotDate]," spread=",string[rw`spread]," pnl=",string[rw`dailyPnl]," cum=",string rw`cumulativePnl;
    srIdx+:1];
-1 "";

/ Black-76 pricing
-1 "--- Black-76 Pricing ---";
callPx:.commodity.black76.price[`call;75f;75f;0.25;0.30;0.05];
putPx:.commodity.black76.price[`put;75f;75f;0.25;0.30;0.05];
callGreeks:.commodity.black76.greeks[`call;75f;75f;0.25;0.30;0.05];
-1 "  ATM call: $",string[callPx],"  put: $",string putPx;
-1 "  delta=",string[callGreeks`delta],"  gamma=",string[callGreeks`gamma],"  vega=",string callGreeks`vega;
-1 "";

/ Kirk spread option
spreadParams:`fwd1`fwd2`strike`expiry`vol1`vol2`correlation`riskFreeRate!(75f;70f;3f;0.25;0.30;0.25;0.8;0.05);
kirkCall:.commodity.spread.kirkPrice[`call;spreadParams];
-1 "--- Kirk Spread Option ---";
-1 "  F1=75 F2=70 K=3: call=$",string kirkCall;
-1 "";

/ Spark spread
-1 "--- Spark Spread Example ---";
sparkSpread:.commodity.electricity.sparkSpread[45f;3.5;7f];
-1 "  power=$45  gas=$3.5  HR=7  spark=$",string sparkSpread;
-1 "";

-1 "=============================================================================";
-1 " Commodity demo completed";
-1 "=============================================================================";
