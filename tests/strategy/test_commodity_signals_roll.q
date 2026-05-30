\l lib/init.q
/ Roll-adjusted returns: hand 2-contract case. Front A rolls to B 5 days before
/ A's expiry; the return on the roll boundary must be B's OWN return (no level
/ jump from A to B). frontContractYM records the held sequence.
contractA:202002; contractB:202003;
expA:2020.02.20; expB:2020.03.20;
dates:2020.02.13 2020.02.14 2020.02.15 2020.02.16 2020.02.17;
pricesA:60 61 62 63 64f;
pricesB:50 50.5 51 51.5 52f;
mkDate:{[i;dates;ca;cb;ea;eb;pa;pb]
    d:dates i;
    ([] asofDate:d,d; tenor:((ea-d)%365f),(eb-d)%365f; price:(pa i),pb i; contractYM:ca,cb; expiry:ea,eb)};
hist:raze mkDate[;dates;contractA;contractB;expA;expB;pricesA;pricesB] each til count dates;

sigCfg:`rollDaysBeforeExpiry`trainEndDate`momentumLookback`kalmanEstCfg!(
    5;2020.02.15;2;`gridSteps`refineRounds`nSweeps!(3;1;1));
out:.strategy.path.commoditySignals[hist;sigCfg];
p:out`path;

.testutil.assertTrue[5=count p;"5 dated path rows"];
/ Held sequence: A held into dates 0,1,2; roll at date2 -> B held into dates 3,4.
.testutil.assertTrue[(202002 202002 202002 202003 202003)~p`frontContractYM;"held sequence rolls A->B at the trigger"];
/ Return on the roll boundary (index 3) is B's own return, not an A->B jump.
.testutil.assertNear[(p`frontReturn)3;(51.5%51f)-1f;1e-12;"roll-boundary return is the new contract's own return"];
/ The pre-roll return (index 2) is A's own return.
.testutil.assertNear[(p`frontReturn)2;(62f%61f)-1f;1e-12;"pre-roll return is the old contract's own return"];
.testutil.assertTrue[0f=(p`frontReturn)0;"first return is zero (no prior)"];

-1 "PASS test_commodity_signals_roll";
