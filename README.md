<img src="xpresso_logo.png" width="300">

# Xpresso: Predicting gene expression levels from genomic sequence

This repository is intended to accompany our publication, primarily to enhance the reproducibility of our results. For more information please refer to:

Agarwal V, Shendure J. [Predicting mRNA abundance directly from genomic sequence using deep convolutional neural networks](https://www.biorxiv.org/content/10.1101/416685v2). 2020. **_Cell Reports_**.

These tools can be used in a variety of organisms and cell types of interest to:

* Perform hyperparameter optimization in the gene expression prediction task (as shown in Fig 1)
* Perform evolutionary analyses on human and mouse organisms, as well as one-to-one orthologs of each (as shown in Fig 2)
* Uncover modes of gene regulation in a cell type of interest that are operating at the transcriptional and post-transcriptional levels (as shown in Fig 3)
* Evaluate model performance for cell type-specifc and cell type-agnostic models (as shown in Fig 4)
* Predict transcriptional activity across a genomic locus (as shown in Fig 5)
* Interpret deep learning models to learn about promoter properties (as shown in Fig 6)

If you find our code or predictions to be helpful for your work, please cite the paper above.


# Dependencies for running entire pipeline:
* Python3 modules: numpy, h5py, pandas, sklearn, keras (>=2.2.4-tf), hyperopt, biopython

* R libraries: LSD, data.table, latticeExtra, Biostrings, rhdf5, ROCR, gplots, mixtools, reshape2, beeswarm, RColorBrewer, zoo, GenomicRanges

* [TensorFlow (>=1.15.0)](https://www.tensorflow.org/install/)

* [DeepExplain](https://github.com/marcoancona/DeepExplain)

* [The MEME Suite](http://meme-suite.org/doc/download.html?man_type=web)

* [UCSC tools](http://hgdownload.soe.ucsc.edu/downloads.html#source_downloads) installation, including bigBedToBed

* [BEDTools](https://github.com/arq5x/bedtools2/releases)

# Instructions for use

For R code to work properly, please copy the contents of .Rprofile in this folder to your local .Rprofile.

Users are advised to read the code closely and modify commented pieces as appropriate to acquire
desired output for your environment. For example, you will need to download all of the additional
R library and Python module dependencies for the code to work. This being said, if you find crucial
files are missing, making the code unusable, or if you identify a major problem in the code, please
raise a Github issue.

In each Figure's folder, change directories to it and please read the file "runme.sh" first as it provides a general overview of relevant commands that were used sequentially to pre-process the data and generate the figures.

Run the following command in the base Xpresso directory to download the associated datapack:

`wget -r -np -nH --reject "index.html*" --cut-dirs 5 https://krishna.gs.washington.edu/content/members/vagar/Xpresso/data/datasets/`

The figures will link to this folder accordingly. Some of the files need to be decompressed, and not all files are provided due to minimize the package size (currently ~11Gb). If you need additional files not provided for the purpose of reproduction, please contact Vikram Agarwal (vagar {at} calicolabs {dot} com).

# Colab

Start training models and generating predictions quickly using the iPython Notebook,
or open it in Google Colab to get up to use a cloud GPU with this link:
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/gist/vagarwal87/bdd33e66fa2c59c41409ca47e7132e61/xpresso.ipynb)

**Note: The Colab generates predictions on a FASTA file of arbitrary DNA sequences without considering mRNA half-life features. To consider half-life features, one must prepare the full test file as shown
in the datapack and Fig1_S2/**
