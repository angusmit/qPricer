/ test_portfolio_monte_carlo_products.q - portfolio with equity + Asian
\l lib/init.q

/ All trades must have the same columns to form a valid table
equityTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`averageType`averagingStyle`observationCount`barrierType`barrierLevel`rebate!(
    1;`AAPL;`equityOption;`european;`call;100f;1f;1000000f;`none;`none;0N;`none;0Nf;0f);
asianCallTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`averageType`averagingStyle`observationCount`barrierType`barrierLevel`rebate!(
    2;`AAPL;`asianOption;`european;`call;100f;1f;1000000f;`arithmetic;`discrete;50;`none;0Nf;0f);
americanAsianTrade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional`averageType`averagingStyle`observationCount`barrierType`barrierLevel`rebate!(
    3;`AAPL;`asianOption;`american;`call;100f;1f;1000000f;`arithmetic;`discrete;50;`none;0Nf;0f);

tradeTable:();
tradeTable:tradeTable,enlist equityTrade;
tradeTable:tradeTable,enlist asianCallTrade;
tradeTable:tradeTable,enlist americanAsianTrade;

mkt:`underlying`spot`riskFreeRate`dividendYield`volatility!(`AAPL;100f;0.05;0f;0.2);
mdl:.model.createBlackScholesModel[];
mcConfig:.montecarlo.defaultMcConfig[];
cfg:.config.defaultPricingConfig[];
cfgWithMC:cfg,enlist[`mcConfig]!enlist mcConfig;

pricingResult:.portfolio.priceTrades[tradeTable;mkt;mdl;cfgWithMC];

/ 1. Three rows
.testutil.assertTrue[3=count pricingResult;"3 pricing rows"];

/ 2. Equity option OK
equityRow:pricingResult 0;
.testutil.assertTrue[equityRow[`status]~`OK;"equity option priced OK"];
.testutil.assertTrue[equityRow[`unitPrice]>0f;"equity price positive"];

/ 3. Asian call OK
asianRow:pricingResult 1;
.testutil.assertTrue[asianRow[`status]~`OK;"Asian call priced OK"];
.testutil.assertTrue[asianRow[`unitPrice]>0f;"Asian price positive"];

/ 4. American Asian should error
americanRow:pricingResult 2;
.testutil.assertTrue[americanRow[`status]~`ERROR;"American Asian produces error"];

/ 5. notionalPrice = unitPrice * notional for OK rows
.testutil.assertNear[equityRow`notionalPrice;equityRow[`unitPrice]*1000000f;1f;"equity notional"];

-1 "PASS test_portfolio_monte_carlo_products";
