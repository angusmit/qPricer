/ test_sabr_smile_generation.q
\l lib/init.q
sabrParams:`alpha`beta`rho`nu!(0.2;0.5;-0.3;0.4);
strikeList:80 90 100 110 120f;

/ Smile table
smileResult:.sabr.smileTable[100f;strikeList;1f;sabrParams];
.testutil.assertTrue[5=count smileResult;"5 smile rows"];
okCount:sum (smileResult`status)=`OK;
.testutil.assertTrue[okCount=5;"all 5 OK"];
allPositive:all (smileResult`impliedVolatility)>0f;
.testutil.assertTrue[allPositive;"all vols positive"];

/ Surface table
expiryList:0.25 0.5 1f;
surfaceResult:.sabr.surfaceTable[100f;strikeList;expiryList;sabrParams];
.testutil.assertTrue[15=count surfaceResult;"15 surface rows (5*3)"];

-1 "PASS test_sabr_smile_generation: smileRows=",string[count smileResult],", surfaceRows=",string count surfaceResult;
