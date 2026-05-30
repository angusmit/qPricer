\l core/init.q
/ Curve shape sanity:
/   (a) every futures price strictly positive,
/   (b) bumping shortFactor0 lifts near tenors more than far tenors (decay exp(-kappa tau)),
/   (c) bumping longFactor0 shifts every tenor by the same factor exp(deltaY).

params:`meanReversionSpeed`shortVolatility`longVolatility`longDrift`correlation!(1.5;0.35;0.15;0f;0.3);
shortFactor0:0f;
longFactor0:log 75f;
tenorVector:0.25 0.5 1 2 5f;

baseCurve:.commodity.schwartz2.futuresCurve[shortFactor0;longFactor0;params;tenorVector];
.testutil.assertTrue[all baseCurve>0f;"curve strictly positive"];

bumpedShortCurve:.commodity.schwartz2.futuresCurve[shortFactor0+1f;longFactor0;params;tenorVector];
shortLogDiffs:log[bumpedShortCurve]-log baseCurve;
nearShortDiff:first shortLogDiffs;
farShortDiff:last shortLogDiffs;
.testutil.assertTrue[nearShortDiff>farShortDiff;"short factor bump fades with maturity"];

kappaVal:params`meanReversionSpeed;
expectedNearShortDiff:exp neg kappaVal*first tenorVector;
expectedFarShortDiff:exp neg kappaVal*last tenorVector;
.testutil.assertNear[nearShortDiff;expectedNearShortDiff;1e-10;"short bump near matches exp(-k tau_near)"];
.testutil.assertNear[farShortDiff;expectedFarShortDiff;1e-10;"short bump far matches exp(-k tau_far)"];

bumpAmount:0.1;
bumpedLongCurve:.commodity.schwartz2.futuresCurve[shortFactor0;longFactor0+bumpAmount;params;tenorVector];
longLogDiffs:log[bumpedLongCurve]-log baseCurve;
.testutil.assertTrue[(max longLogDiffs)-min longLogDiffs<1e-10;"long factor bump is parallel"];
.testutil.assertNear[first longLogDiffs;bumpAmount;1e-10;"long factor bump magnitude matches deltaY"];

-1 "PASS test_schwartz2_futures_curve_shape: baseCurve=",-3!baseCurve;
