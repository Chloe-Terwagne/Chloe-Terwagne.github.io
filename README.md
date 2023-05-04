# VarEffectViz
VarEffectViz is a visualizer that presents variant effect interpretation for the BRCA1 gene using multiple sources of data, including allele count in UKBiobank, the CADD computational score, AlphaFold protein structure prediction and the results of saturation genome editing experiments.



The SGE function score in the vizualisation is the SGE function score is the normalized negative function score obtained using min-max normalization on function score from Findlay et al., 2018

Cumulative score is computed as follow:

$cumulative\_{score} = \frac{1}{3} \left(\frac{1}{AC} + SGE_{score} + \frac{cadd_{score}}{50}\right)$

where $AC$ is allele count, $SGE_{score}$ is the normalized negative function score obtained using min-max normalization on function score from Findlay et al., 2018, and cadd_score is the CADD score divided by the maximum cadd score in this region which is 50.

Note that 1/3 is a constant factor used to normalize the sum of the three scores.
