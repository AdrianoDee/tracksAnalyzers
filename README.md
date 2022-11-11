# tracksAnalyzers

Usage:

``` cmsRun trackAnalyzerAOD.py f=AOD_INPUT_FILE T=THREADS e=NUMBEROFEVENTS```

``` cmsRun trackAnalyzerMINIAOD.py f=MINIAOD_INPUT_FILE T=THREADS e=NUMBEROFEVENTS```

produces a flat root tree with all the tracks above 0.6 GeV. The  `cov_matrix.ipynb` pour them into a pandas daframe ad makes some plots for the covariance.
