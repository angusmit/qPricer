/ limitcheck.q - limit check comparison utilities (v0.26)

.limitcheck.absoluteDifference:{[leftPrice;rightPrice] abs leftPrice-rightPrice};

.limitcheck.relativeDifference:{[leftPrice;rightPrice]
    avgPrice:0.5*abs[leftPrice]+abs rightPrice;
    if[avgPrice<1e-12; :0f];
    (abs leftPrice-rightPrice)%avgPrice
 };

.limitcheck.withinTolerance:{[leftPrice;rightPrice;absoluteTolerance;relativeTolerance]
    absDiff:.limitcheck.absoluteDifference[leftPrice;rightPrice];
    relDiff:.limitcheck.relativeDifference[leftPrice;rightPrice];
    (absDiff<=absoluteTolerance) or relDiff<=relativeTolerance
 };

.limitcheck.resultRow:{[checkName;leftModel;rightModel;leftPrice;rightPrice;absoluteTolerance;relativeTolerance]
    absDiff:.limitcheck.absoluteDifference[leftPrice;rightPrice];
    relDiff:.limitcheck.relativeDifference[leftPrice;rightPrice];
    passedVal:.limitcheck.withinTolerance[leftPrice;rightPrice;absoluteTolerance;relativeTolerance];
    `checkName`leftModel`rightModel`leftPrice`rightPrice`absoluteDifference`relativeDifference`absoluteTolerance`relativeTolerance`passed`status`errorMessage!(
        checkName;leftModel;rightModel;leftPrice;rightPrice;absDiff;relDiff;absoluteTolerance;relativeTolerance;passedVal;`OK;"")
 };

.limitcheck.errorRow:{[checkName;leftModel;rightModel;errMsg]
    `checkName`leftModel`rightModel`leftPrice`rightPrice`absoluteDifference`relativeDifference`absoluteTolerance`relativeTolerance`passed`status`errorMessage!(
        checkName;leftModel;rightModel;0Nf;0Nf;0Nf;0Nf;0Nf;0Nf;0b;`ERROR;errMsg)
 };

.limitcheck.summary:{[checkResultTable]
    if[0=count checkResultTable;
        :`checkCount`passedCount`failedCount`maxAbsoluteDifference`maxRelativeDifference`status!(0;0;0;0Nf;0Nf;`OK)];
    passedCol:checkResultTable`passed;
    passedCount:sum passedCol;
    totalCount:count checkResultTable;
    failedCount:totalCount-passedCount;
    absDiffs:checkResultTable`absoluteDifference;
    relDiffs:checkResultTable`relativeDifference;
    validAbs:absDiffs where not null absDiffs;
    validRel:relDiffs where not null relDiffs;
    `checkCount`passedCount`failedCount`maxAbsoluteDifference`maxRelativeDifference`status!(
        totalCount;passedCount;failedCount;
        $[0<count validAbs;max validAbs;0Nf];
        $[0<count validRel;max validRel;0Nf];
        `OK)
 };
