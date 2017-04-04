# Dissociation curve of the H2O dimer, and CDFT example

In this experiment, we use the fragment approximation for the dissociation curve of a H2O dimer, and we generalize the
results for a Constrained DFT calculation of a H2O - OH- dimer.

## H2O dimer

The goal of this exercice is to test the fragment approach of BigDFT and to compare the results against the one obtained using the cubic and the linear versions of BigDFT. This is the last exercise of the BigDFT fragment tutorial. The studied problem is the dissociation curve of a H2O dimer, the key parameter being the distance between the two water molecules. The influence of different input parameters are studied.


## H2O - HO-

For this second experiment, the system under study has the same number of electrons as before, but is lacking a hydrogen atom. Cubic simulations have been done as well as fragment approach simulations (no full linear simulations yet).
The goal of this part of the notebook is to illustrate how to use constrained DFT in the fragment approach with BigDFT. Two initial setups were used: one using OH as a template, the other using OH- as the template. Constrained DFT can be applied on top of that, in order to force the added electron to remain in the OH fragment.
