\l core/init.q
/ ============================================================================
/ test_directional_template_faithful.q - the byte-identity proof. The directionalSignal
/ template DELEGATES to the unchanged signal+backtest path, so its per-day PnL must EQUAL
/ the existing backtest's per-day PnL exactly. Synthetic equity path (gammaScalp), no HDB.
/ Research OS R6. This is the test that keeps the abstraction honest - it can't silently
/ change the directional numbers.
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};

pathCfg:`spot0`drift`volatility`steps`stepYears`riskFreeRate`dividendYield`seed!(100f;0f;0.30;11;1f%252f;0.05;0f;42);
pathTbl:.strategy.path.fromSynthetic pathCfg;
trade:`tradeId`underlying`productType`exerciseStyle`optionType`strike`expiry`notional!(`T_GS;`X;`equityOption;`european;`call;100f;0.25;1f);
bsModel:.model.createBlackScholesModel[];
fdmCfg:`method`numberOfSpotSteps`numberOfTimeSteps`minimumSpot`maximumSpot`interpolationMethod`returnFullGrid`stabilityCheck!(`crankNicolson;100;200;0f;300f;`linear;1b;1b);
stratCfg:@[.strategy.defaultConfig `gammaScalp;(`rebalanceMode;`rebalanceInterval;`stepYears);:;(`interval;1;1f%252f)];

/ the EXISTING direct path
direct:select date:stepDate, pnl:stepPnl from
    (.strategy.runAndSummarize[`gammaScalp;trade;pathTbl;bsModel;fdmCfg;stratCfg])`result where status=`OK;

/ the SAME run via the directionalSignal template (a faithful wrapper)
tmpl:.template.run[`directionalSignal;
    `strategy`trade`path`model`fdmCfg`stratCfg!(`gammaScalp;trade;pathTbl;bsModel;fdmCfg;stratCfg)];

chk[direct~tmpl`pnl; "the directional template's per-day PnL must be BYTE-IDENTICAL to the direct backtest"];
chk[(exec sum pnl from direct)~exec sum pnl from tmpl`pnl; "and therefore the total PnL is identical"];
chk[(tmpl[`meta;`shape])~`directionalSignal; "the template meta must record its shape"];
chk[(tmpl[`meta;`strategy])~`gammaScalp; "the template meta must record the wrapped strategy"];

/ a second strategy (short-gamma config) - still byte-identical
shortCfg:@[stratCfg;`optionSide;:;`short];
direct2:select date:stepDate, pnl:stepPnl from
    (.strategy.runAndSummarize[`gammaScalp;trade;pathTbl;bsModel;fdmCfg;shortCfg])`result where status=`OK;
tmpl2:.template.run[`directionalSignal; `strategy`trade`path`model`fdmCfg`stratCfg!(`gammaScalp;trade;pathTbl;bsModel;fdmCfg;shortCfg)];
chk[direct2~tmpl2`pnl; "a different config is still wrapped byte-identically"];

-1 "test_directional_template_faithful: PASS";
