/ correlation.q - correlation matrix construction, validation, Cholesky (v0.17)

/ --- Public ---

.correlation.validateCorrelationTable:{[correlationTable;symbolList]
    numSymbols:count symbolList;
    if[numSymbols<2; '"Need at least 2 symbols for correlation"];
    requiredPairCount:(numSymbols*(numSymbols-1))%2;
    if[not requiredPairCount=count correlationTable;
        '"Expected ",string[requiredPairCount]," correlation pairs for ",string[numSymbols]," symbols, got ",string count correlationTable];
    corrValues:correlationTable`correlation;
    if[any corrValues>1f; '"Correlation values must be <= 1"];
    if[any corrValues< -1f; '"Correlation values must be >= -1"];
 };

.correlation.toMatrix:{[correlationTable;symbolList]
    numSymbols:count symbolList;
    matrixData:numSymbols cut (numSymbols*numSymbols)#1f;
    pairIdx:0;
    while[pairIdx<count correlationTable;
        pairRow:correlationTable pairIdx;
        sym1Idx:symbolList?pairRow`sym1;
        sym2Idx:symbolList?pairRow`sym2;
        if[(sym1Idx>=numSymbols) or sym2Idx>=numSymbols;
            '"Symbol not found in symbolList: ",(string pairRow`sym1)," or ",string pairRow`sym2];
        corrVal:pairRow`correlation;
        matrixData[sym1Idx;sym2Idx]:corrVal;
        matrixData[sym2Idx;sym1Idx]:corrVal;
        pairIdx+:1];
    matrixData
 };

.correlation.isSymmetric:{[matrixData;toleranceVal]
    transposed:flip matrixData;
    maxDiff:max max each abs matrixData-transposed;
    maxDiff<=toleranceVal
 };

.correlation.isPositiveSemiDefinite:{[matrixData;toleranceVal]
    / Try Cholesky — if it succeeds, matrix is PD
    choleskyResult:@[.correlation.__cholesky;matrixData;{`FAIL}];
    not choleskyResult~`FAIL
 };

.correlation.choleskyOrNearPsd:{[matrixData]
    .correlation.__cholesky matrixData
 };

.correlation.nearestPsd:{[matrixData]
    / Simple approach: add epsilon to diagonal until PD
    dimSize:count matrixData;
    adjustedMatrix:matrixData;
    epsilon:1e-10;
    attemptIdx:0;
    while[attemptIdx<20;
        choleskyResult:@[.correlation.__cholesky;adjustedMatrix;{`FAIL}];
        if[not choleskyResult~`FAIL; :choleskyResult];
        idx:0;
        while[idx<dimSize; adjustedMatrix[idx;idx]:adjustedMatrix[idx;idx]+epsilon; idx+:1];
        epsilon:epsilon*10f;
        attemptIdx+:1];
    '"Failed to find nearest PSD after 20 attempts"
 };

/ --- Internal: Cholesky decomposition ---

.correlation.__cholesky:{[matrixData]
    dimSize:count matrixData;
    choleskyL:dimSize cut (dimSize*dimSize)#0f;
    rowIdx:0;
    while[rowIdx<dimSize;
        colIdx:0;
        while[colIdx<=rowIdx;
            if[colIdx=rowIdx;
                sumSquared:0f;
                kIdx:0;
                while[kIdx<colIdx;
                    sumSquared+:choleskyL[rowIdx;kIdx]*choleskyL[rowIdx;kIdx];
                    kIdx+:1];
                diagVal:matrixData[rowIdx;colIdx]-sumSquared;
                if[diagVal<=0f; '"Cholesky failed: matrix is not positive definite"];
                choleskyL[rowIdx;colIdx]:sqrt diagVal];
            if[not colIdx=rowIdx;
                dotProduct:0f;
                kIdx:0;
                while[kIdx<colIdx;
                    dotProduct+:choleskyL[rowIdx;kIdx]*choleskyL[colIdx;kIdx];
                    kIdx+:1];
                choleskyL[rowIdx;colIdx]:(matrixData[rowIdx;colIdx]-dotProduct)%choleskyL[colIdx;colIdx]];
            colIdx+:1];
        rowIdx+:1];
    choleskyL
 };
