/ =============================================================================
/ Example: Validate real market data against qPricer FDM
/ =============================================================================
/ Reads options data from CSV and compares market prices to FDM model prices
/ Run from project root: q examples/validate_market_data.q
/ =============================================================================

\l core/init.q

/ ----------------------------
/ Configuration
/ ----------------------------

/ Path to CSV file (edit this for your data)
csvPath:hsym `$"C:/Users/angus/Downloads/aapl-options-exp-2026-05-26-weekly-near-the-money-stacked-05-25-2026.csv";

/ Market snapshot date and values
valDate:2026.05.25;
spotPrice:308.82;
riskFreeRate:0.045;
dividendYield:0f;

/ FDM solver configuration
fdmConfig:`method`numberOfSpotSteps`numberOfTimeSteps!(
    `explicit;
    200;
    2000
 );

/ ----------------------------
/ Helper Functions
/ ----------------------------

/ Parse percentage string (e.g., "23.45%" -> 0.2345)
pctToFloat:{[x]
    "F"$ssr/[string x;("+";"%";"unch";"N/A");("";"";"0";"0")]
 };

/ Parse number string (e.g., "unch" -> 0)
numToFloat:{[x]
    "F"$ssr/[string x;("unch";"N/A");("0";"0")]
 };

/ Price a single option using qPricer engine
priceOption:{[trade;marketData]
    model:.model.createBlackScholesModel[];
    config:.config.createFiniteDifferenceConfig[fdmConfig];
    result:.engine.priceOption[trade;marketData;model;config];
    result`unitPrice
 };

/ Create trade record from market data row
mkTrade:{[row;valDate]
    optTypeVal:$[row`optType~`call;`call;`put];
    `tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(
        .tradeIdCounter+:1;
        row`symbol;
        `equityOption;
        `european;
        optTypeVal;
        row`strike;
        (row`expDate-valDate)%365f;
        1f
    )
 };

/ Create market data record from row
mkMarketData:{[row;valDate;spot;rate;qy]
    `underlying`spot`riskFreeRate`dividendYield`volatility!(
        row`symbol;
        spot;
        rate;
        qy;
        row`impliedVol
    )
 };

/ ----------------------------
/ Main: Load & Validate
/ ----------------------------

-1 "=== qPricer Market Data Validation ===";
-1 "Loading: ",string csvPath;

.tradeIdCounter:0;

/ Read CSV
raw:-1 _ read0 csvPath;
tbl:("FSFFFFSSIISSFSS"; enlist ",") 0: raw;

/ Column mapping
tbl:`strike`moneyness`bid`mid`ask`lastPx`chg`pctChg`volume`openInt`oiChange`impliedVol`delta`optType`lastTrade xcol tbl;

/ Clean data
tbl:update
    impliedVol:(pctToFloat each impliedVol)%100,
    pctChg:(pctToFloat each pctChg)%100,
    chg:numToFloat each chg,
    symbol:`AAPL,
    expDate:2026.05.26,
    optType:`$optType
    from tbl;

/ Filter valid rows
valid:select from tbl where bid>0, ask>bid, mid>0, impliedVol>0, volume>0, openInt>0;

-1 "Total rows: ",(string count tbl);
-1 "Valid rows: ",(string count valid);
-1 "";

/ Price each option
.trades:();
results:();


{[rowIdx]
    row:valid[rowIdx];

    trade:mkTrade[row;valDate];
    mdata:mkMarketData[row;valDate;spotPrice;riskFreeRate;dividendYield];

    modelPrice:priceOption[trade;mdata];

    .trades,:enlist trade;
    results,:enlist `symbol`optType`strike`expDate`marketMid`modelPrice`impliedVol`error!(
        row`symbol;
        row`optType;
        row`strike;
        row`expDate;
        row`mid;
        modelPrice;
        row`impliedVol;
        modelPrice - row`mid
    );
 } each til count valid;

results:5!results;

/ ----------------------------
/ Output Results
/ ----------------------------

-1 "=== Validation Results ===";
show select symbol,optType,strike,marketMid,modelPrice,impliedVol,error
    from results;

-1 "";
-1 "Summary by Option Type:";
show select avgPrice:avg modelPrice, avgError:avg error, maxError:max abs error
    by optType from results;

-1 "";
-1 "Worst model errors:";
show 5#`error xdesc select symbol,optType,strike,marketMid,modelPrice,error
    from results;

/ ----------------------------
/ Save Results
/ ----------------------------

outPath:hsym `$"C:/Users/angus/Downloads/qprcer_validation_results";
outPath set results;
-1 "Saved to: ",string outPath;

