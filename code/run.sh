# 物种名字？
# 纠正符号，与黎良的人一致
# 富集从5到5
# motifRank从4到5

###############################################################################
# Usage: bash run.sh
# The script used for running other scripts all in sequence.
###############################################################################

source ~/.bashrc
conda activate r43

# matlab -nodisplay -nosplash -nodesktop < ./SVD/robustSVD.m

Rscript ./motifRank.R

Rscript ./PCC.R

Rscript ./projection.R

# Rscript ./enrichmentAnalysis/gene_enrichment_analysis.R

# Rscript ./enrichmentAnalysis/format_enrichment_result.R
