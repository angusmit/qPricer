/ test_scenario_risk.q — validate scenario risk report
\l core/init.q

trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1000000f);
marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(
    `AAPL;100f;0.05;0f;0.2);
model:.model.createBlackScholesModel[];
config:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

riskReport:.risk.generateScenarioReport[trade;marketData;model;config];

/ 1. Check row count
if[not 9=count riskReport; '"FAIL test_scenario_risk: expected 9 rows, got ",string count riskReport];

/ 2. Check all scenario names present
expectedScenarios:`base`spotUp1Pct`spotDown1Pct`spotUp5Pct`spotDown5Pct`volatilityUp1Point`volatilityDown1Point`rateUp25Bp`rateDown25Bp;
if[not all expectedScenarios in riskReport`scenario; '"FAIL test_scenario_risk: missing scenarios"];

/ 3. Base PnL = 0
baseRow:riskReport where riskReport[`scenario]=`base;
if[not 0f=baseRow[`unitPnL] 0; '"FAIL test_scenario_risk: base unitPnL not zero"];
if[not 0f=baseRow[`notionalPnL] 0; '"FAIL test_scenario_risk: base notionalPnL not zero"];

/ 4. Directional checks for call
baseUnitPrice:baseRow[`unitPrice] 0;
spotUpPrice:(riskReport where riskReport[`scenario]=`spotUp1Pct)[`unitPrice] 0;
spotDnPrice:(riskReport where riskReport[`scenario]=`spotDown1Pct)[`unitPrice] 0;
volUpPrice:(riskReport where riskReport[`scenario]=`volatilityUp1Point)[`unitPrice] 0;
volDnPrice:(riskReport where riskReport[`scenario]=`volatilityDown1Point)[`unitPrice] 0;
rateUpPrice:(riskReport where riskReport[`scenario]=`rateUp25Bp)[`unitPrice] 0;
rateDnPrice:(riskReport where riskReport[`scenario]=`rateDown25Bp)[`unitPrice] 0;

if[not spotUpPrice > baseUnitPrice; '"FAIL: spotUp should increase call price"];
if[not spotDnPrice < baseUnitPrice; '"FAIL: spotDown should decrease call price"];
if[not volUpPrice > baseUnitPrice; '"FAIL: volUp should increase call price"];
if[not volDnPrice < baseUnitPrice; '"FAIL: volDown should decrease call price"];
if[not rateUpPrice > baseUnitPrice; '"FAIL: rateUp should increase call price"];
if[not rateDnPrice < baseUnitPrice; '"FAIL: rateDown should decrease call price"];

/ 5. notionalPnL consistency
pnlTolerance:1e-6;
unitPnLs:riskReport`unitPnL;
notionalPnLs:riskReport`notionalPnL;
expectedNotionalPnLs:unitPnLs * trade`notional;
maxPnLDiff:max abs notionalPnLs - expectedNotionalPnLs;
if[maxPnLDiff > pnlTolerance; '"FAIL: notionalPnL inconsistent with unitPnL * notional"];

spotUpPnL:(riskReport where riskReport[`scenario]=`spotUp1Pct)[`unitPnL] 0;
volUpPnL:(riskReport where riskReport[`scenario]=`volatilityUp1Point)[`unitPnL] 0;
rateUpPnL:(riskReport where riskReport[`scenario]=`rateUp25Bp)[`unitPnL] 0;

-1 "PASS test_scenario_risk: rows=",string[count riskReport],", basePrice=",string[baseUnitPrice],", spotUpPnL=",string[spotUpPnL],", volUpPnL=",string[volUpPnL],", rateUpPnL=",string rateUpPnL;
