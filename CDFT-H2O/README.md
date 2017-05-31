# Dissociation curves of the water dimer and the H2O + HO- dimer: using fragments and Constrained DFT

This set of four notebooks aims is a presents how to use the fragment approach and Constrained DFT with BigDFT. 

In the first two notebooks, the system under study is the water dimer, and the quantity of interest is the dissociation curve (*i.e.* the energy of the dimer as a function of the distance between both molecules). The use of the fragment approach to tackle this problem is presented and compared to the standard cubic and linear calculations.

In the last two notebooks, the system is very similar, except that an hydrogen atom is removed from a water molecule and an extra charge is added. The goal is still to compute the dissociation curve. The fragment approach (and its comparison to the standard cubic and linear calculations) is still considered, but the main point consists in applying constrained DFT to constrain the extra charge to stay in the vicinity of a choosen fragment.

It is advised to run the notebooks in the order given below.

## H2O dimer

### H2O\_dimer.ipynb

The first notebook presents how to use python to launch calculations and how to read and treat the BigDFT logfiles, taking the example of the fragment tutorial on the BigDFT website. 

The fragment approach is used to compute the dissociation curve of the water dimer. We present how to solve the last exercise, which consists in comparing the dissociation curves obtained using the fragment approach with an inscreasing number of support functions to the standard cubic and linear BigDFT calculations.


### H2O\_dimer\_tests.ipynb

This second notebook reviews the influence of other input parameters that may have an impact on the BigDFT calculations, and therefore on the dissociation curves. The main parameters concern the grid extension, the localization radii, the confining potential and finally the possibility to optimize the templates support functions.


## H2O - HO- dimer

### H2O\_HO-.ipynb

This notebook is very similar to the first notebook, the main difference being the system under study, which is the water - hydroxide dimer descibed in the introduction. The influence of the number of support functions per atom is presented.

In contrast with the previous calculations, there are two templates to be used, one for the water molecule, the other for the hydroxide molecule. While templates previously calculated for the water dimer study can be reused for the water molecule, new templates have to be computed for the hydroxide molecule. We propose to use either the HO- template (with an explicit charge added to the template calculation) or the HO template (without the explicit charge). Note that whichever the template used, an extra charge has to be added to the full system.

The dissociation curves obtained using the fragment approach are compared to the reference cubic calculation. 


### H2O\_HO-\_CDFT.ipynb

This notebook is almost the same as the previous one, except that this time the extra charge added to the full system is constrained to remain in the vicinity of the HO fragment. Both HO- and HO templates are used, and the number of support functions per atom is increased.
