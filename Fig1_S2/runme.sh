
#run for ~1000 iterations on a fast GPU (can take ~1-2 days to run)
python deep_motif_learning_hyperopt.py train trials/tpe_1K_10epochs_optimized_0to20K.hyperopt pM10Kb_1KTest/

Rscript Fig1B_S2A.R FILE1 FILE2

for x in {1..10}; do { python deep_motif_learning_hyperopt.py test trials/tpe_1K_10epochs_optimized_0to20K.hyperopt pM10Kb_1KTest/ >trial_$x.txt; } done &

Rscript Fig1C.R FILE1