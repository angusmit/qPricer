/ pnl.q - PnL explain via Taylor expansion using Greeks
/ actualPnL = pv1 - pv0
/ explainedPnL = deltaPnL + gammaPnL + vegaPnL + rhoPnL + thetaPnL
/ unexplainedPnL = actualPnL - explainedPnL

/ --- Public ---

.pnl.explainTrade:{[tradeRow;marketDataBook0;marketDataBook1;configDict]
    pricingModel:configDict`model;
    fdmConfig:configDict`fdmConfig;
    timeStepYears:configDict`timeStepYears;
    bookLabel:configDict`bookName;
    / Get market data for trade from both books
    mktData0:.marketbook.getMarketDataForTrade[marketDataBook0;tradeRow];
    mktData1:.marketbook.getMarketDataForTrade[marketDataBook1;tradeRow];
    / Price at both dates
    priceResult0:.engine.priceOption[tradeRow;mktData0;pricingModel;fdmConfig];
    priceResult1:.engine.priceOption[tradeRow;mktData1;pricingModel;fdmConfig];
    pv0Val:priceResult0`unitPrice;
    pv1Val:priceResult1`unitPrice;
    actualPnLVal:pv1Val-pv0Val;
    / Greeks at t0
    greeksResult:.greeks.calculateGreeks[tradeRow;mktData0;pricingModel;fdmConfig];
    deltaVal:greeksResult[`delta]0;
    gammaVal:greeksResult[`gamma]0;
    thetaVal:greeksResult[`theta]0;
    vegaVal:greeksResult[`vega]0;
    rhoVal:greeksResult[`rho]0;
    / Market data changes
    dSpot:mktData1[`spot]-mktData0`spot;
    dVol:mktData1[`volatility]-mktData0`volatility;
    dRate:mktData1[`riskFreeRate]-mktData0`riskFreeRate;
    / PnL components
    deltaPnLVal:deltaVal*dSpot;
    gammaPnLVal:0.5*gammaVal*dSpot*dSpot;
    vegaPnLVal:vegaVal*dVol;
    rhoPnLVal:rhoVal*dRate;
    thetaPnLVal:thetaVal*timeStepYears;
    explainedPnLVal:deltaPnLVal+gammaPnLVal+vegaPnLVal+rhoPnLVal+thetaPnLVal;
    unexplainedPnLVal:actualPnLVal-explainedPnLVal;
    `tradeId`bookName`underlyingSym`pv0`pv1`actualPnL`deltaPnL`gammaPnL`vegaPnL`rhoPnL`thetaPnL`explainedPnL`unexplainedPnL`status`errorMessage!(
        tradeRow`tradeId;bookLabel;tradeRow`underlying;
        pv0Val;pv1Val;actualPnLVal;
        deltaPnLVal;gammaPnLVal;vegaPnLVal;rhoPnLVal;thetaPnLVal;
        explainedPnLVal;unexplainedPnLVal;`OK;"")
 };

.pnl.explainPortfolio:{[tradeTable;marketDataBook0;marketDataBook1;configDict]
    numTrades:count tradeTable;
    resultList:();
    loopIdx:0;
    while[loopIdx<numTrades;
        currentTrade:tradeTable loopIdx;
        singleResult:.[.pnl.explainTrade;(currentTrade;marketDataBook0;marketDataBook1;configDict);{x}];
        if[10h=type singleResult;
            singleResult:`tradeId`bookName`underlyingSym`pv0`pv1`actualPnL`deltaPnL`gammaPnL`vegaPnL`rhoPnL`thetaPnL`explainedPnL`unexplainedPnL`status`errorMessage!(
                currentTrade`tradeId;configDict`bookName;currentTrade`underlying;
                0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;`ERROR;singleResult)];
        resultList:resultList,enlist singleResult;
        loopIdx+:1];
    resultList
 };

.pnl.aggregateExplain:{[explainResult]
    okRows:explainResult where explainResult[`status]=`OK;
    if[0=count okRows;
        :`totalPv0`totalPv1`totalActualPnL`totalDeltaPnL`totalGammaPnL`totalVegaPnL`totalRhoPnL`totalThetaPnL`totalExplainedPnL`totalUnexplainedPnL!(
            0f;0f;0f;0f;0f;0f;0f;0f;0f;0f)];
    `totalPv0`totalPv1`totalActualPnL`totalDeltaPnL`totalGammaPnL`totalVegaPnL`totalRhoPnL`totalThetaPnL`totalExplainedPnL`totalUnexplainedPnL!(
        sum okRows`pv0;sum okRows`pv1;sum okRows`actualPnL;
        sum okRows`deltaPnL;sum okRows`gammaPnL;sum okRows`vegaPnL;
        sum okRows`rhoPnL;sum okRows`thetaPnL;sum okRows`explainedPnL;sum okRows`unexplainedPnL)
 };
