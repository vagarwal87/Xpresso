# precomputed h5 files for human (pM10Kb_cv) and mouse (pM10Kb_Mouse_cv) are provided only to save space
# but all can be generated as below:
python setup_training_files.py --cv Roadmap_FantomAnnotations.InputData.pM10Kb.txt.gz pM10Kb_cv
python setup_training_files.py --cv Roadmap_FantomAnnotations.InputData.pM10Kb.K562expr.txt.gz pM10Kb_1KTest_K562expr_cv
python setup_training_files.py --cv Roadmap_FantomAnnotations.InputData.pM10Kb.GM12878expr.txt.gz pM10Kb_1KTest_GM12878expr_cv
python setup_training_files.py --cv Mouse_FantomAnnotations.InputData.pM10Kb.txt.gz pM10Kb_Mouse_cv
python setup_training_files.py --cv Mouse_FantomAnnotations.InputData.pM10Kb.mESC.txt.gz pM10Kb_1KTest_mESCexpr_cv

# RUN ON GPU USING FOLDERS ABOVE, TAKES MANY HOURS TO RUN ON GPU
for y in {1..10}; do { for x in {0..9}; do { python Xpresso.py --fold $y --trial $x test tpe_1K_10epochs_optimized_0to20K.hyperopt pM10Kb_cv/ >pM10Kb_cv/fold$y\_trial$x.txt; } done } done &
for y in {1..10}; do { for x in {0..9}; do { python Xpresso.py --fold $y --trial $x test tpe_1K_10epochs_optimized_0to20K.hyperopt pM10Kb_1KTest_K562expr_cv/ >pM10Kb_1KTest_K562expr_cv/fold$y\_trial$x.txt; } done } done &
for y in {1..10}; do { for x in {0..9}; do { python Xpresso.py --fold $y --trial $x test tpe_1K_10epochs_optimized_0to20K.hyperopt pM10Kb_1KTest_GM12878expr_cv/ >pM10Kb_1KTest_GM12878expr_cv/fold$y\_trial$x.txt; } done } done &
for y in {1..10}; do { for x in {0..9}; do { python Xpresso.py --fold $y --trial $x test tpe_1K_10epochs_optimized_0to20K.hyperopt pM10Kb_Mouse_cv/ >pM10Kb_Mouse_cv/fold$y\_trial$x.txt; } done } done &
for y in {1..10}; do { for x in {0..9}; do { python Xpresso.py --fold $y --trial $x test tpe_1K_10epochs_optimized_0to20K.hyperopt pM10Kb_1KTest_mESCexpr_cv/ >pM10Kb_1KTest_mESCexpr_cv/fold$y\_trial$x.txt; } done } done &

#MERGED RESULTS INTO ALL-CROSSVALIDATED PREDICTIONS
Rscript integrate_cv_results.R cross_valid pM10Kb_cv all_crossvalidated_predictions.txt
Rscript integrate_cv_results.R cross_valid_K562 pM10Kb_1KTest_K562expr_cv all_crossvalidated_predictions_K562.txt
Rscript integrate_cv_results.R cross_valid_GM12878 pM10Kb_1KTest_GM12878expr_cv all_crossvalidated_predictions_GM12878.txt
Rscript integrate_cv_results.R cross_valid_mouse pM10Kb_Mouse_cv all_crossvalidated_predictions_mouse.txt
Rscript integrate_cv_results.R cross_valid_mESC pM10Kb_1KTest_mESCexpr_cv all_crossvalidated_predictions_mESC.txt

mkdir diHMM
cd diHMM
wget http://bcb.dfci.harvard.edu/~gcyuan/data/diHMM/diHMM_Annotations.zip
unzip diHMM_Annotations.zip

cd K562
bedtools intersect -wo -a K562_nD30_nB30_domainLevelStatesColor.bed -b ../../Homo_sapiens.hg19.90.chosenTranscript.geneBoundaries.gtf.gz >K562_overlapping_genes.bed
grep -P 'D7|D8|D9|D23' K562_overlapping_genes.bed | cut -f 18 | cut -b 1-15 | sort | uniq >H3K27me3_silenced.txt
grep -P 'D10|D11|D12|D13' K562_overlapping_genes.bed | cut -f 18 | cut -b 1-15 | sort | uniq >superenhancer.txt
cd ../GM12878/
bedtools intersect -wo -a GM12878_nD30_nB30_domainLevelStatesColor.bed -b ../../Homo_sapiens.hg19.90.chosenTranscript.geneBoundaries.gtf.gz >GM12878_overlapping_genes.bed
grep -P 'D7|D8|D9|D23' GM12878_overlapping_genes.bed | cut -f 18 | cut -b 1-15 | sort | uniq >H3K27me3_silenced.txt
grep -P 'D10|D11|D12|D13' GM12878_overlapping_genes.bed | cut -f 18 | cut -b 1-15 | sort | uniq >superenhancer.txt
cd ../..

Rscript Fig3ABCDEF_S3ABC.R

wget http://www.targetscan.org/mmu_71/mmu_71_data_download/Summary_Counts.default_predictions.txt.zip
unzip Summary_Counts.default_predictions.txt.zip
gzip Summary_Counts.default_predictions.txt.gz
Rscript Fig3GH.R
