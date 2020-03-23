# Deep Explain CV computes the following methods:
# 'deeplift', 'grad*input', 'saliency', 'elrp', 'intgrad'
# 'grad*input' and 'intgrad' are presented in the paper for the reasons described
# thus only these two have stored precomputed results to save space in download

#human (using best model from 10 trials on each of 10 folds of data)
python deep_explain_cv.py cv_human/01trainepoch.10-0.427.h5 pM10Kb_cv/ 1
python deep_explain_cv.py cv_human/62trainepoch.13-0.421.h5 pM10Kb_cv/ 2
python deep_explain_cv.py cv_human/63trainepoch.10-0.428.h5 pM10Kb_cv/ 3
python deep_explain_cv.py cv_human/24trainepoch.09-0.427.h5 pM10Kb_cv/ 4
python deep_explain_cv.py cv_human/25trainepoch.16-0.422.h5 pM10Kb_cv/ 5
python deep_explain_cv.py cv_human/76trainepoch.11-0.414.h5 pM10Kb_cv/ 6
python deep_explain_cv.py cv_human/17trainepoch.09-0.414.h5 pM10Kb_cv/ 7
python deep_explain_cv.py cv_human/08trainepoch.11-0.418.h5 pM10Kb_cv/ 8
python deep_explain_cv.py cv_human/79trainepoch.14-0.436.h5 pM10Kb_cv/ 9
python deep_explain_cv.py cv_human/110trainepoch.08-0.425.h5 pM10Kb_cv/ 10

#Fig6A and FigS6A-B
Rscript FigS6.R pM10Kb_cv/ human #resulting plots in {human/mouse}{gradinput/intgrad}.pdf

for x in {1..10}; do { Rscript best_positions.R pM10Kb_cv/gradinput.$x.txt.gz $x human; } done

cd motif_analysis/
#generate permuted set
for x in {1..10}; do { cut -f 1 negbestpos1mer$x\human.txt | shuf | paste - <(cut -f 2- bestpos1mer$x\human.txt) >negbestpos1mer$x\human.txt; } done
for x in {1..10}; do { ./extract_kmers.pl human <bestpos1mer$x\human.txt >bestpos1mer$x\human.fa; } done
for x in {1..10}; do { ./extract_kmers.pl human <negbestpos1mer$x\human.txt >negbestpos1mer$x\human.fa; } done
cat bestpos1mer*human.fa >bestpos1mer_human_all.fa
cat negbestpos1mer*human.fa >negbestpos1mer_human_all.fa
fasta-subsample bestpos1mer_human_all.fa 400000 >subsampled_bestpos1mer_human_all_400K.fa
fasta-subsample negbestpos1mer_human_all.fa 400000 >subsampled_negbestpos1mer_human_all_400K.fa
#Fig6B human
dreme -oc human_all_1mer_400K/ -p subsampled_bestpos1mer_human_all_400K.fa -n subsampled_negbestpos1mer_human_all_400K.fa -mink 2 -maxk 10
cd ..

#mouse (using best model from 10 trials on each of 10 folds of data)
python deep_explain_cv.py cv_mouse/71trainepoch.07-0.3200.h5 pM10Kb_Mouse_cv/ 1
python deep_explain_cv.py cv_mouse/02trainepoch.07-0.3186.h5 pM10Kb_Mouse_cv/ 2
python deep_explain_cv.py cv_mouse/63trainepoch.06-0.3173.h5 pM10Kb_Mouse_cv/ 3
python deep_explain_cv.py cv_mouse/64trainepoch.09-0.3194.h5 pM10Kb_Mouse_cv/ 4
python deep_explain_cv.py cv_mouse/45trainepoch.13-0.3113.h5 pM10Kb_Mouse_cv/ 5
python deep_explain_cv.py cv_mouse/96trainepoch.07-0.3134.h5 pM10Kb_Mouse_cv/ 6
python deep_explain_cv.py cv_mouse/77trainepoch.07-0.3223.h5 pM10Kb_Mouse_cv/ 7
python deep_explain_cv.py cv_mouse/88trainepoch.06-0.3243.h5 pM10Kb_Mouse_cv/ 8
python deep_explain_cv.py cv_mouse/79trainepoch.07-0.3171.h5 pM10Kb_Mouse_cv/ 9
python deep_explain_cv.py cv_mouse/610trainepoch.10-0.3200.h5 pM10Kb_Mouse_cv/ 10

#FigS6C-D
Rscript FigS6.R pM10Kb_Mouse_cv/ mouse

for x in {1..10}; do { Rscript best_positions.R pM10Kb_Mouse_cv/gradinput.$x.txt.gz $x mouse; } done

cd motif_analysis/
#generate permuted set
for x in {1..10}; do { cut -f 1 negbestpos1mer$x\mouse.txt | shuf | paste - <(cut -f 2- bestpos1mer$x\mouse.txt) >negbestpos1mer$x\mouse.txt; } done
for x in {1..10}; do { ./extract_kmers.pl mouse <bestpos1mer$x\mouse.txt >bestpos1mer$x\mouse.fa; } done
for x in {1..10}; do { ./extract_kmers.pl mouse <negbestpos1mer$x\mouse.txt >negbestpos1mer$x\mouse.fa; } done
cat bestpos1mer*mouse.fa >bestpos1mer_mouse_all.fa
cat negbestpos1mer*mouse.fa >negbestpos1mer_mouse_all.fa
fasta-subsample bestpos1mer_mouse_all.fa 400000 >subsampled_bestpos1mer_mouse_all_400K.fa
fasta-subsample negbestpos1mer_mouse_all.fa 400000 >subsampled_negbestpos1mer_mouse_all_400K.fa
#Fig6B mouse
dreme -oc mouse_all_1mer_400K/ -p subsampled_bestpos1mer_mouse_all_400K.fa -n subsampled_negbestpos1mer_mouse_all_400K.fa -mink 2 -maxk 10
cd ..

# for human (Fig6C and S7)
Rscript Fig6C_S7_S8.R Roadmap_FantomAnnotations.InputData.pM10Kb.txt.gz FigS7.pdf
# for mouse (FigS8)
Rscript Fig6C_S7_S8.R Mouse_FantomAnnotations.InputData.pM10Kb.txt.gz FigS8.pdf
