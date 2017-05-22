# Profiling examples of BigDFT timings with futile python module

For a given production run of BigDFT, a file named "time[-]<run_name>.yaml" is created close to the logfile (or in the "data[-]<run_name>/" directory). Such a file might be used to investigate the main features of the run and can be used for comparisons with other files for strong and weak scaling analysis.
We provide in this directory few examples on how such files can be used to infer the code behaviour on various machines.

We suggest to start by having a look to the file "PerformanceInvestigation.ipynb" first as it describes the basis for the analysis.
