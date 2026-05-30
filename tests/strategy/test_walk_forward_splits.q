\l lib/init.q
/ Walk-forward splitter: known-answer ranges (non-overlapping test windows, no
/ leak: trainEnd < testEnd) and INDEX-STABILITY (appending future data adds new
/ splits without changing existing ones - the cross-split causality guarantee).
expSplits:.strategy.commodityBT.__splits[`expanding;30;10;5;8];
.testutil.assertTrue[4=count expSplits;"expanding: 4 splits fit in 30 dates"];
.testutil.assertTrue[(0 0 0 0)~expSplits[;`trainStartIdx];"expanding train always starts at 0"];
.testutil.assertTrue[(9 14 19 24)~expSplits[;`trainEndIdx];"expanding trainEnd grows by testSpan"];
.testutil.assertTrue[(14 19 24 29)~expSplits[;`testEndIdx];"expanding testEnd = trainEnd + testSpan"];
/ No leak: each split's trainEnd strictly precedes its testEnd.
.testutil.assertTrue[all expSplits[;`trainEndIdx]<expSplits[;`testEndIdx];"trainEnd < testEnd (no leak)"];

rollSplits:.strategy.commodityBT.__splits[`rolling;30;10;5;8];
.testutil.assertTrue[(0 5 10 15)~rollSplits[;`trainStartIdx];"rolling train start slides by testSpan"];
.testutil.assertTrue[all 10=1+rollSplits[;`trainEndIdx]-rollSplits[;`trainStartIdx];"rolling train width is constant (10)"];

/ Index stability: the first splits are identical when the history grows.
bigSplits:.strategy.commodityBT.__splits[`expanding;50;10;5;8];
.testutil.assertTrue[expSplits~(count expSplits)#bigSplits;"existing splits unchanged when data is appended"];

/ A test window that runs past the data is not emitted.
.testutil.assertTrue[0=count .strategy.commodityBT.__splits[`expanding;12;10;5;8];"no split when test window exceeds the data"];

-1 "PASS test_walk_forward_splits";
