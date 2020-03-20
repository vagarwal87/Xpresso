Rscript Fig2EFG.R #required to run before Fig 2A-D

### RAN FROM DATASETS DIRECTORY ON GPU, OUTPUT HERE FOR CONVENIENCE
# mkdir subsampled
# for y in {0..9}; do { for x in `seq 2000 2000 16000`; do { python deep_motif_learning_hyperopt.py -c $x test trials/tpe_1K_10epochs_optimized_0to20K.hyperopt pM10Kb_1KTest/ >subsampled/train$y\_$x.txt; } done } done &

Rscript Fig2A.R
Rscript Fig2BC.R

#ORTHOLOG cmd here ortholog_results/
Rscript Fig2D.R
