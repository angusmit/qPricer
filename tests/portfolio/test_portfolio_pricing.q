/ test_portfolio_pricing.q - portfolio batch pricing
\l lib/init.q

tradeTable:([]
    tradeId:1 2 3 4;
    underlying:`AAPL`AAPL`AAPL`AAPL;
    productType:`equityOption`equityOption`equityOption`equityOption;
    exerciseStyle:`european`european`american`european;
    optionType:`call`put`put`call;
    strike:100 100 100 100f;
    expiry:1 1 1 1f;
    notional:1000000 1000000 1000000 1000000f;
    barrierType:`none`none`none`upAndOut;
    barrierLevel:0N 0N 0N 130f;
    rebate:0 0 0 0f);

marketData:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
model:.model.createBlackScholesModel[];
config:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(
    `explicit;200;2000;0f;300f;`linear;1b;1b);

portfolioPriceTable:.portfolio.priceTrades[tradeTable;marketData;model;config];

/ 1. Correct row count
if[not 4=count portfolioPriceTable; '"FAIL: expected 4 rows"];

/ 2. All status OK
statusValues:portfolioPriceTable`status;
if[not all statusValues=`OK; '"FAIL: not all status OK"];

/ 3. All prices positive
priceValues:portfolioPriceTable`unitPrice;
if[any priceValues<=0f; '"FAIL: non-positive unitPrice"];

/ 4. All notional prices positive
notionalValues:portfolioPriceTable`notionalPrice;
if[any notionalValues<=0f; '"FAIL: non-positive notionalPrice"];

/ 5. European call ~ 10.455
euroCallPrice:priceValues 0;
if[(abs euroCallPrice-10.45496)>0.01; '"FAIL: European call price unexpected"];

/ 6. European put ~ 5.578
euroPutPrice:priceValues 1;
if[(abs euroPutPrice-5.577839)>0.01; '"FAIL: European put price unexpected"];

/ 7. American put > European put
americanPutPrice:priceValues 2;
if[not americanPutPrice>euroPutPrice; '"FAIL: American put should exceed European put"];

/ 8. Barrier call < vanilla call
barrierCallPrice:priceValues 3;
if[not barrierCallPrice<euroCallPrice; '"FAIL: barrier call should be less than vanilla call"];

totalNotional:sum notionalValues;
-1 "PASS test_portfolio_pricing: rows=4, totalNotionalPrice=",string totalNotional;
