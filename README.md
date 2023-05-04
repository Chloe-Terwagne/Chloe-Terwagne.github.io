# VarEffectViz
VarEffectViz is a visualizer that presents variant effect interpretation for the BRCA1 gene using multiple sources of data, including allele count in UKBiobank, the CADD computational score, AlphaFold protein structure prediction and the results of saturation genome editing experiments.

Cumulative score is computed as follow:
cumulative_score = (1/3) * (1/AC + minmax_neg_func_score + cadd_score/50)

where AC is allele count, minmax_neg_func_score is the normalized negative function score obtained using min-max normalization, and cadd_score is the CADD score divided by the maximum cadd score in this region which is 50.

Note that 1/3 is a constant factor used to normalize the sum of the three scores.
