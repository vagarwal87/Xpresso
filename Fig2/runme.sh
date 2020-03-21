### Ran from datasets/ directory on GPU
mkdir subsampled_10fold
for x in `seq 2000 2000 16000`; do { echo $x; python subsample.py $x pM10Kb_1KTest subsampled_10fold/$x; } done
for z in `seq 2000 2000 16000`; do { for y in {0..9}; do { for x in {0..9}; do { python Xpresso.py --fold $y --trial $x train tpe_1K_10epochs_optimized_0to20K.hyperopt subsampled_10fold/$z/ >subsampled_10fold/$z/$y\_trial$x.txt; } done } done } done

Rscript Fig2A.R
Rscript Fig2BC.R

#Acquired 1 to 1 ortholog predictions from Ensembl BioMart
grep one2one human2mouse_orthologs.txt >human2mouse_one2one_orthologs.txt
python setup_training_files.py -t 1000 -v 1000 --orthologs Mouse_FantomAnnotations.InputData.pM10Kb.txt.gz Roadmap_FantomAnnotations.InputData.pM10Kb.txt.gz pM10Kb_1KTest_one2oneOrthologs

mkdir ortholog_results
mkdir ortholog_results/train_human_test_human
mkdir ortholog_results/train_mouse_test_human
mkdir ortholog_results/train_human_test_mouse
mkdir ortholog_results/train_mouse_test_mouse
for x in {1..10}; do { python Xpresso.py train tpe_1K_10epochs_optimized_0to20K.hyperopt pM10Kb_1KTest_one2oneOrthologs/ >ortholog_results/train_human_test_human/trial_$x.txt; } done &
for x in {1..10}; do { python Xpresso.py train tpe_1K_10epochs_optimized_0to20K.hyperopt pM10Kb_1KTest_one2oneOrthologs/ >ortholog_results/train_human_test_mouse/trial_$x.txt; } done &
for x in {1..10}; do { python Xpresso.py train tpe_1K_10epochs_optimized_0to20K.hyperopt pM10Kb_1KTest_one2oneOrthologs/ >ortholog_results/train_mouse_test_human/trial_$x.txt; } done &
for x in {1..10}; do { python Xpresso.py train tpe_1K_10epochs_optimized_0to20K.hyperopt pM10Kb_1KTest_one2oneOrthologs/ >ortholog_results/train_mouse_test_mouse/trial_$x.txt; } done &

#Stored results from runs in ortholog_results/
Rscript Fig2D.R
Rscript Fig2EFG.R
