# Copyright © 2024, Jianhui Shi & Lei M. Li. Academy of Mathematics and Systems Science, Chinese Academy of Sciences, Beijing 100190, China

###############################################################################
# This script is used for running other scripts all in sequence.
###############################################################################

matlab -nodisplay -nosplash -nodesktop < ./SVD/robustSVD.m

Rscript ./motifRank.R

Rscript ./PCC.R

Rscript ./projection.R

Rscript ./enrichmentAnalysis/gene_enrichment_analysis.R

Rscript ./enrichmentAnalysis/format_enrichment_result.R
