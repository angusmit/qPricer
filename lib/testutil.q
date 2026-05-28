/ testutil.q - reusable test assertion helpers

.testutil.assertTrue:{[condition;testMessage]
    if[not condition; '"ASSERT FAILED: ",testMessage];
 };

.testutil.assertEqual:{[actualVal;expectedVal;testMessage]
    if[not actualVal~expectedVal;
        '"ASSERT FAILED: ",testMessage," | expected: ",(-3!expectedVal)," | actual: ",(-3!actualVal)];
 };

.testutil.assertNear:{[actualVal;expectedVal;toleranceVal;testMessage]
    if[(abs actualVal-expectedVal)>toleranceVal;
        '"ASSERT FAILED: ",testMessage," | expected~",string[expectedVal]," actual=",string[actualVal]," tol=",string toleranceVal];
 };

.testutil.assertTableColumns:{[tbl;requiredColumns;testMessage]
    tableColumnNames:cols tbl;
    missingColumns:requiredColumns where not requiredColumns in tableColumnNames;
    if[0<count missingColumns;
        '"ASSERT FAILED: ",testMessage," | missing: ",", " sv string missingColumns];
 };

.testutil.expectError:{[testName;functionToRun]
    functionResult:@[{x[];`NO_ERROR};functionToRun;{`ERROR}];
    if[functionResult~`NO_ERROR; '"Expected error but got success: ",testName];
    -1 "  PASS ",testName;
 };
