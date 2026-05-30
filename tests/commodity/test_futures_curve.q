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

.commodity.futures.validateFuturesTable futuresTable;
-1 "  PASS validate";

futuresCurve:.commodity.futures.buildCurveByDate futuresTable;
.testutil.assertTrue[(count futuresCurve)>0;"curve has rows"];
.testutil.assertTrue[`tenorRank in key futuresCurve 0;"has tenorRank"];

/ Front contract
frontRow:.commodity.futures.frontContract[futuresCurve;2024.01.02];
.testutil.assertTrue[frontRow[`tenorRank]=1;"front is tenor 1"];
.testutil.assertNear[frontRow`settlementPrice;72.5;0.01;"front price"];

/ Nth contract
secondRow:.commodity.futures.nthContract[futuresCurve;2024.01.02;2];
.testutil.assertTrue[secondRow[`tenorRank]=2;"second tenor"];

/ Curve slope
slopeVal:.commodity.futures.curveSlope[futuresCurve;2024.01.02;1;3];
.testutil.assertNear[slopeVal;1.3;0.01;"slope = far - near"];

/ Invalid table
badResult:@[.commodity.futures.validateFuturesTable;();{`ERROR}];
.testutil.assertTrue[badResult~`ERROR;"empty table fails"];

-1 "PASS test_futures_curve";
