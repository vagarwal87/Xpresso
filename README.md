<img src="xpresso_logo.png" width="200">

# Xpresso: Predicting gene expression levels from genomic sequence

This repository is intended to accompany our manuscript and enhance the reproducibility of our results. For more information please refer to:

Agarwal V, Shendure J. [Predicting mRNA abundance directly from genomic sequence using deep convolutional neural networks](https://www.biorxiv.org/content/10.1101/416685v2). _**Cell Reports**_. (2020).

These tools can be used in a variety of organisms and cell types of interest to:

* Perform hyperparameter optimization in the gene expression prediction task (as shown in Fig 1)
* Perform evolutionary analyses on human and mouse organisms, as well as one-to-one orthologs of each (as shown in Fig 2)
* Uncover modes of gene regulation in a cell type of interest that are operating at the transcriptional and post-transcriptional levels (as shown in Fig 3)
* Evaluate model performance for cell type-specifc and cell type-agnostic models (as shown in Fig 4)
* Predict transcriptional activity across a genomic locus (as shown in Fig 5)
* Interpret deep learning models (as shown in Fig 6)

If you find our code or predictions to be helpful for your work, please cite the paper above.


# Dependencies for running entire pipeline:
* Python modules: numpy, h5py, pandas, sklearn, keras (tested on v2.0.8), hyperopt, tensorflow, biopython

* R libraries: LSD, data.table, latticeExtra, Biostrings, rhdf5, ROCR, gplots

* [TensorFlow (tested on v1.3.0)] (https://www.tensorflow.org/install/)

* [The MEME Suite] (http://meme-suite.org/doc/download.html?man_type=web)

* [UCSC tools](http://hgdownload.soe.ucsc.edu/downloads.html#source_downloads) installation, including bigBedToBed

* [BEDTools] (https://github.com/arq5x/bedtools2/releases)

# Instructions for use

For R code to work properly, please copy the contents of .Rprofile in this folder to your local .Rprofile.

Users are advised to read the code closely and modify commented pieces as appropriate to acquire
desired output for your environment. For example, you will need to download all of the additional
R library and Python module dependencies for the code to work. This being said, if you find crucial
files are missing, making the code unusable, or if you identify a major problem in the code, please
raise a Github issue.

In each Figure's folder, change directories to it and please read the file "runme.sh" first as it provides a general overview of relevant commands that were used sequentially to pre-process the data and generate the figures.

Run the following command in the base Xpresso directory to download the associated datapack and data
preparation scripts:

`wget -r -np -nH --reject "index.html*" --cut-dirs 5 \
 https://krishna.gs.washington.edu/content/members/vagar/Xpresso/data/datasets/`

# Additional notes
