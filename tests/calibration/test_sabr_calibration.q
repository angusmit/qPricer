/ test_sabr_calibration.q
\l core/init.q

/ Generate synthetic market smile from known params
trueParams:`alpha`beta`rho`nu!(0.2;0.5;-0.3;0.4);
strikeList:80 90 100 110 120f;
syntheticSmile:.sabr.smileTable[100f;strikeList;1f;trueParams];

/ Build market smile table
smileTable:();
smileIdx:0;
while[smileIdx<count syntheticSmile;
    smileRow:syntheticSmile smileIdx;
    smileTable:smileTable,enlist `strike`expiry`marketImpliedVolatility`weight!(
        smileRow`strike;smileRow`expiry;smileRow`impliedVolatility;1f);
    smileIdx+:1];

/ Grid including true params
paramGrid:`alphaList`betaList`rhoList`nuList!(
    0.15 0.2 0.25;
    enlist 0.5;
    -0.5 -0.3 0.0;
    0.2 0.4 0.6);

sabrGrid:.sabr.calibrationGrid paramGrid;
calibTable:.sabr.calibrateSmile[smileTable;100f;sabrGrid;()!()];
.testutil.assertTrue[(count calibTable)>0;"calibration has rows"];

bestRow:.sabr.bestCalibration calibTable;
.testutil.assertTrue[bestRow[`rmse]<0.01;"best RMSE small"];
.testutil.assertNear[bestRow`alpha;0.2;0.06;"best alpha near true"];
.testutil.assertNear[bestRow`rho;-0.3;0.21;"best rho near true"];

-1 "PASS test_sabr_calibration: bestRmse=",string[bestRow`rmse],", bestAlpha=",string[bestRow`alpha],", bestRho=",string bestRow`rho;
