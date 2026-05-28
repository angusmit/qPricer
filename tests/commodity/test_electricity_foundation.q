\l lib/init.q
/ Synthetic hub prices
hubTable:();
hIdx:0;
basePrices:45 48 120 46 44 50 47;
hubDates:2024.01.01+til 7;
while[hIdx<7;
    hubTable:hubTable,enlist `snapshotDate`hub`region`price`volume!(hubDates hIdx;`PJM;`east;basePrices hIdx;10000f);
    hIdx+:1];

.commodity.electricity.validateHubPriceTable hubTable;

/ Seasonality features
withSeason:.commodity.electricity.seasonalityFeatures hubTable;
.testutil.assertTrue[`monthOfYear in key withSeason 0;"has monthOfYear"];
.testutil.assertTrue[`dayOfWeek in key withSeason 0;"has dayOfWeek"];

/ Spike detection
withSpikes:.commodity.electricity.spikeFlags[hubTable;enlist[`nStd]!enlist 2f];
.testutil.assertTrue[`spikeFlag in key withSpikes 0;"has spikeFlag"];
/ 120 should be a spike
spikeRow:withSpikes 2;
.testutil.assertTrue[spikeRow`spikeFlag;"120 is spike"];

/ Monthly average
monthAvg:.commodity.electricity.monthlyAveragePrice hubTable;
.testutil.assertTrue[(count monthAvg)>0;"monthly avg rows"];

/ Spark spread
gasTable:();
gIdx:0;
gasPrices:3.5 3.6 3.4 3.5 3.3 3.7 3.5;
while[gIdx<7;
    gasTable:gasTable,enlist `snapshotDate`hub`region`price`volume!(hubDates gIdx;`HenryHub;`us;gasPrices gIdx;5000f);
    gIdx+:1];

sparkResult:.commodity.electricity.sparkSpreadTable[hubTable;gasTable;7f];
.testutil.assertTrue[`sparkSpread in key sparkResult 0;"has sparkSpread"];
/ 45 - 7*3.5 = 45 - 24.5 = 20.5
.testutil.assertNear[(sparkResult 0)`sparkSpread;20.5;0.01;"spark spread day1"];

-1 "PASS test_electricity_foundation";
