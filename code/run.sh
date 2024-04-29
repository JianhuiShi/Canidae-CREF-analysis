###############################################################################
# This script is used for running other scripts all in sequence.
###############################################################################

matlab -nodisplay -nosplash -nodesktop < ./SVD/robustSVD.m

Rscript ./motifRank.R

Rscript ./PCC.R

Rscript ./projection.R

Rscript ./enrichmentAnalysis/gene_enrichment_analysis.R

Rscript ./enrichmentAnalysis/format_enrichment_result.R
