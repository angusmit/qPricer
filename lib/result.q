/ result.q - standardised result/error handling
/ Uses `OK/`ERROR to match existing portfolio/pnl/batch convention.

.result.okRow:{[dataDict]
    dataDict,`status`errorMessage!(`OK;"")
 };

.result.errorRow:{[dataDict;errorMsg]
    dataDict,`status`errorMessage!(`ERROR;errorMsg)
 };

.result.isOk:{[resultTable]
    resultTable[`status]=`OK
 };

.result.isError:{[resultTable]
    not resultTable[`status]=`OK
 };

.result.errorSummary:{[resultTable]
    resultTable where not resultTable[`status]=`OK
 };

.result.assertColumns:{[resultTable;requiredColumns]
    tableColumnNames:cols resultTable;
    missingColumns:requiredColumns where not requiredColumns in tableColumnNames;
    if[0<count missingColumns;
        '"Missing required columns: ",", " sv string missingColumns];
    1b
 };
