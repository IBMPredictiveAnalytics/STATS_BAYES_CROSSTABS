﻿* Encoding: UTF-8.
file handle data /name="c:/spss24/samples/english".

get file = "data/employee data.sav".
dataset name emp.
compute randomwt = rv.uniform(.01, 1).

stats bayes crosstabs rows=educ columns=jobcat /options priorconcentration=100.



STATS BAYES CROSSTABS ROWS=educ COLUMNS=jobcat SAMPLETYPE=POISSON 
/OPTIONS DISPLAYTABLE=YES POSTERIOR=YES PRIORCONCENTRATION=1
/SAVE WORKSPACE=CLEAR.

weight by randomwt.

STATS BAYES CROSSTABS ROWS=educ COLUMNS=jobcat SAMPLETYPE=POISSON 
/OPTIONS DISPLAYTABLE=YES POSTERIOR=YES PRIORCONCENTRATION=1
/SAVE WORKSPACE=CLEAR.

DATASET ACTIVATE emp.
STATS BAYES CROSSTABS ROWS=educ COLUMNS=jobcat SAMPLETYPE=POISSON 
/OPTIONS DISPLAYTABLE=YES POSTERIOR=YES PRIORCONCENTRATION=.1
/SAVE WORKSPACE=CLEAR.
weight off.
