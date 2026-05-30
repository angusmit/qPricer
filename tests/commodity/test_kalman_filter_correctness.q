\l lib/init.q
/ Kalman filter correctness on a trivial KNOWN linear-Gaussian state space (a
/ scalar local level: alpha_t = alpha_{t-1} + eta, y_t = alpha_t + eps). The
/ generic .commodity.kalman.filterBundles must match an INDEPENDENT hand-coded
/ scalar Kalman recursion, in both the filtered mean and the log-likelihood.
processVar:0.5; measVar:0.2; ys:1.0 0.5 2.0;

/ 1x1 matrix bundles for the local level.
mkBundle:{[yObs;q;r]
    `Tm`c`W`Z`d`V`y!((enlist enlist 1f);(enlist 0f);(enlist enlist q);(enlist enlist 1f);(enlist 0f);(enlist enlist r);(enlist yObs))};
bundles:mkBundle[;processVar;measVar] each ys;
a0:(enlist 0f); P0:(enlist enlist 1f);
fb:.commodity.kalman.filterBundles[bundles;a0;P0];
filterMeans:(fb`states)[;`a][;0];

/ Independent scalar Kalman recursion (closed-form local level).
twoPi:2f*acos neg 1f;
handStep:{[st;yObs]
    a:st 0; pVar:st 1; ll:st 2;
    aPred:a; pPred:pVar+0.5;
    innov:yObs-aPred; sVar:pPred+0.2;
    kGain:pPred%sVar; aNew:aPred+kGain*innov; pNew:pPred-kGain*pPred;
    llNew:ll-0.5*(log[2f*acos neg 1f]+log[sVar]+innov*innov%sVar);
    (aNew;pNew;llNew)};
handStates:handStep\[(0f;1f;0f);ys];
handMeans:handStates[;0];
handLoglik:(last handStates)2;

.testutil.assertTrue[1e-10>abs (fb`loglik)-handLoglik;"filter log-likelihood matches hand recursion"];
.testutil.assertTrue[1e-10>max abs filterMeans-handMeans;"filtered means match hand recursion"];
.testutil.assertTrue[3=count filterMeans;"one filtered state per observation"];

-1 "PASS test_kalman_filter_correctness: loglik=",string fb`loglik;
