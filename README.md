# partialDistanceCorrelation

This matlab package provides several ways to perform partial distance correlation. The scripts in the main directory (mainly pdc.m) perform *partial* bias corrected distance correlation between two variables of interest and any number of variables that the relationship should be conditioned on. This is achieved by applying the recursive formula for partial correlation, but to bias corrected distance correlation coefficients. The recursive formula is implemented to handle any nth order partial correlation. P-values are generated in the same way as for bias corrected distance correlation (chi-square distribution). In addition, a linear partial correlation can be performed using the same recursive formula (i.e., for comparison). This can also be compared to matlab's internal partialcorr (also output by pdc.m), which shows very small (i.e., rounding, precision) deviations with higher order partial correlations (3+). 

You can now use pdcPerm.m to perform a permutation based test of the significance of the bias corrected partial distance correlation (pdc). The y variable is permuted and the p-value is computed as one-tail (i.e., proportion of times permuted pdc coefficient exceeds actual pdc coefficient). In other words, the null is whether the pdc coefficient is higher than chance level pdc. 

You can also use pdcPerm even if you don't want to do a permutation analysis. I recommend this as pdcPerm will clean up NaNs for you and check variables.

Note, pdcPerm uses a different method for computing partial correlation, but will produce the same output as dcor in python (see: https://dcor.readthedocs.io/en/latest/index.html)

For more information, see: 
Székely, G. J., & Rizzo, M. L. (2014). Partial distance correlation with methods for dissimilarities. The Annals of Statistics, 42(6), 2382-2412.
Also see these slides: https://stat.wisc.edu/wp-content/uploads/sites/870/2020/03/SzekelyGabor.pdf

In ./permDCwithPartialCorrelation non-bias corrected partial distance correlation is computed and a permutation analysis is used to obtain p-values. partialdistcorr.m can compute the distance correlation after *linear* regression of control variables from two variables of interest. As such, this package does not consider nonlinear relationships between control variables and variables of interest. 

Alex Teghipco // alex.teghipco@sc.edu
