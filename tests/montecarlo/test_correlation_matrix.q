/ test_correlation_matrix.q - correlation matrix construction and validation
\l core/init.q

.test.expectError:{[testName;fn]
    testResult:@[{x[];`NO_ERROR};fn;{`ERROR}];
    if[testResult~`NO_ERROR; '"Expected error but got success: ",testName];
    -1 "  PASS ",testName;
 };

symbolList:`AAPL`MSFT`NVDA;
corrTable:([] sym1:`AAPL`AAPL`MSFT; sym2:`MSFT`NVDA`NVDA; correlation:0.5 0.3 0.4);

/ 1. Valid table builds matrix
.correlation.validateCorrelationTable[corrTable;symbolList];
corrMatrix:.correlation.toMatrix[corrTable;symbolList];
.testutil.assertTrue[3=count corrMatrix;"3x3 matrix"];

/ 2. Diagonal is 1
.testutil.assertNear[corrMatrix[0;0];1f;0.001;"diagonal 0 is 1"];
.testutil.assertNear[corrMatrix[1;1];1f;0.001;"diagonal 1 is 1"];
.testutil.assertNear[corrMatrix[2;2];1f;0.001;"diagonal 2 is 1"];

/ 3. Symmetric
.testutil.assertTrue[.correlation.isSymmetric[corrMatrix;0.001];"matrix is symmetric"];

/ 4. PSD
.testutil.assertTrue[.correlation.isPositiveSemiDefinite[corrMatrix;0.001];"matrix is PSD"];

/ 5. Cholesky succeeds
choleskyL:.correlation.__cholesky corrMatrix;
.testutil.assertTrue[3=count choleskyL;"Cholesky returns 3x3"];

/ 6. Invalid correlation > 1
.test.expectError["correlation > 1";{
    badTable:([] sym1:enlist `AAPL; sym2:enlist `MSFT; correlation:enlist 1.5);
    .correlation.validateCorrelationTable[badTable;`AAPL`MSFT]}];

/ 7. Non-PSD matrix fails Cholesky
.test.expectError["non-PSD matrix";{
    badMatrix:(1 0.99 0.99;0.99 1 -0.99;0.99 -0.99 1f);
    .correlation.__cholesky badMatrix}];

-1 "PASS test_correlation_matrix";
