\l core/init.q
/ realCurveCalendarRoll drives the existing commodityCalendar over a real-shaped
/ curve history. Structural checks: the curve bundle is built correctly, the
/ driver runs, and the regime breakdown of rollPnl sums to the result total
/ (accounting itself is inherited from commodityCalendar's own tests).
contracts:202002 202003 202004 202005;
expiries:2020.02.20 2020.03.20 2020.04.20 2020.05.20;
dates:2020.01.13+7*til 10;
mkRows:{[d;contracts;expiries]
    alive:where expiries>=d;
    tens:(`float$(expiries alive)-d)%365f;
    base:60f-3f*tens;
    ([] asofDate:(count alive)#d; tenor:tens; price:base; contractYM:contracts alive; expiry:expiries alive)};
hist:raze mkRows[;contracts;expiries] each dates;

/ Bundle builder produces a commodityCalendar-compatible bundle.
bundle:.strategy.path.curveBundleFromHistory hist;
.testutil.assertTrue[all `frontPath`curveSnapshots`tenors`evolutionParams`frontLevels in key bundle;"bundle has commodityCalendar keys"];
.testutil.assertTrue[10=count bundle`frontPath;"frontPath one row per date"];
.testutil.assertNear[first bundle`frontLevels;first hist`price;1e-9;"front level = nearest-tenor price"];

out:.strategy.realCurveCalendarRoll[hist;`nearTenor`farTenor`rollTriggerTenor`stepYears!(0.10;0.25;0.20;7f%365f)];
res:out`result;
.testutil.assertTrue[10=count res;"driver produced one row per date"];
.testutil.assertTrue[all (res`status)=`OK;"all steps OK"];
/ Regime breakdown of rollPnl reconciles to the total over OK rows.
okRes:res where (res`status)=`OK;
.testutil.assertNear[sum (out`rollByRegime)`rollPnlTotal;sum okRes`rollPnl;1e-9;"regime rollPnl sums to total"];
.testutil.assertTrue[all (out`rollByRegime)[`regime] in `backwardation`contango;"regimes classified"];

-1 "PASS test_strategy_real_curve_calendar_roll";
