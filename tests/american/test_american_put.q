/ test_american_put.q — validate American put pricing
/ Key checks:
/   1. American put prices without error
/   2. American put >= European put (early exercise premium >= 0)
/   3. Early exercise premium is positive for ATM put with r > 0
/   4. Deep ITM American put approaches intrinsic value
\l core/init.q

/ --- Standard parameters ---
europeanPutTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    200;`TEST;`equityOption;`european;`put;100f;1f;1f);

americanPutTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
    201;`TEST;`equityOption;`american;`put;100f;1f;1f);

marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(
    `TEST;100f;0.05;0f;0.2);

model:.model.createBlackScholesModel[];
config:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

/ --- Test 1: American put prices without error ---
americanResult:.engine.priceOption[americanPutTrade;marketData;model;config];
americanUnitPrice:americanResult`unitPrice;
if[null americanUnitPrice; '"FAIL: American put returned null price"];

/ --- Test 2: European put for comparison ---
europeanResult:.engine.priceOption[europeanPutTrade;marketData;model;config];
europeanUnitPrice:europeanResult`unitPrice;

/ --- Test 3: American put >= European put ---
earlyExercisePremium:americanUnitPrice - europeanUnitPrice;
if[earlyExercisePremium < -0.001;
    '"FAIL: American put (",string[americanUnitPrice],") < European put (",string[europeanUnitPrice],")"];

/ --- Test 4: Early exercise premium > 0 for ATM put with r > 0 ---
if[not earlyExercisePremium > 0.0;
    '"FAIL: Expected positive early exercise premium, got ",string earlyExercisePremium];

/ --- Test 5: Deep ITM American put approaches intrinsic ---
deepItmMkt:@[marketData;`spot;:;50f];
deepItmResult:.engine.priceOption[americanPutTrade;deepItmMkt;model;config];
deepItmPrice:deepItmResult`unitPrice;
deepItmIntrinsic:100f - 50f;
if[deepItmPrice < deepItmIntrinsic - 0.5;
    '"FAIL: Deep ITM American put (",string[deepItmPrice],") below intrinsic (",string[deepItmIntrinsic],")"];

/ --- Test 6: exerciseStyle in metadata ---
americanGridResult:.engine.priceOptionWithGrid[americanPutTrade;marketData;model;config];
solverMeta:(americanGridResult`solverResult)`metadata;
if[not solverMeta[`exerciseStyle]~`american;
    '"FAIL: metadata exerciseStyle not american"];

-1 "PASS test_american_put: american=",string[americanUnitPrice],", european=",string[europeanUnitPrice],", earlyExercisePremium=",string earlyExercisePremium;
