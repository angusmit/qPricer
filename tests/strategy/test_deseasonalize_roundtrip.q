\l lib/init.q
/ Deseasonalization round-trip + known-answer: dividing each contract's price by
/ exp(monthFactor[deliveryMonth]) then multiplying back recovers the original
/ (round-trip), and a purely-seasonal curve deseasonalizes to flat (the seasonal
/ slope is removed). Independent of the strategy signals.
seasonal:0.15 0.10 0.0 -0.08 -0.12 -0.10 -0.05 -0.03 0.0 0.02 0.06 0.12;
ymList:raze {[y] (y*100)+1+til 12} each 2021 2022;
expiryOf:{[ym] "D"$(string[ym div 100]),".",(-2#"0",string[ym mod 100]),".15"};
expiries:expiryOf each ymList;
mkRow:{[d;ymList;expiries;seasonal]
    alive:where expiries>=d;
    months:(ymList[alive] mod 100)-1;
    ([] asofDate:(count alive)#d; tenor:(`float$(expiries alive)-d)%365f;
        price:3.0*exp seasonal months; contractYM:ymList alive; expiry:expiries alive)};
hist:raze mkRow[;ymList;expiries;seasonal] each 2021.01.05+7*til 30;

/ Round-trip with arbitrary factors.
mf:0.10 0.05 -0.02 0.0 -0.05 -0.10 -0.08 -0.03 0.0 0.02 0.04 0.09;
deseason:.strategy.path.__deseasonalizeCurve[hist;mf];
reseason:update price:price*exp mf .strategy.path.__monthIdxFromYM contractYM from deseason;
.testutil.assertTrue[1e-12>max abs (reseason`price)-hist`price;"deseasonalize then reseasonalize == original"];
.testutil.assertTrue[not (deseason`price)~hist`price;"deseasonalization actually changes prices"];

/ Known-answer: fit on the whole history, deseasonalize, the near-vs-far log slope
/ collapses (a purely-seasonal curve has no genuine carry).
fitted:.strategy.path.__fitSeasonalTrain[hist;last hist`asofDate];
.testutil.assertTrue[12=count fitted;"12 fitted monthly factors"];
desFit:.strategy.path.__deseasonalizeCurve[hist;fitted];
avgAbsSlope:{[ch]
    dts:asc distinct ch`asofDate;
    s:{[ch;d] sub:`tenor xasc select tenor,price from ch where asofDate=d;
        $[2>count sub;0n;((log sub[`price]0)-log sub[`price]1)%(sub[`tenor]1)-sub[`tenor]0]}[ch;] each dts;
    avg abs 0f^s};
.testutil.assertTrue[(avgAbsSlope desFit)<0.1*avgAbsSlope hist;"deseasonalized curve slope collapses vs raw seasonal slope"];

-1 "PASS test_deseasonalize_roundtrip";
