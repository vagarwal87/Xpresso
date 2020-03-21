import sys, os, h5py
import numpy.random as npr
import numpy as np
from optparse import OptionParser
import pandas as pd
from sklearn import preprocessing
from sklearn.model_selection import KFold

def main():
    usage = 'usage: %prog [options] <data_file> <out_dir>'
    parser = OptionParser(usage)
    parser.add_option('-t', dest='testCount', default=1000, type='int', help='Number of test examples: [Default: %default]')
    parser.add_option('-v', dest='validCount', default=1000, type='int', help='Number of validation examples: [Default: %default]')
    parser.add_option('--cv', dest='crossVal', default=False, action='store_true', help='Generate samples for 10-fold cross-validated predictions? [Default: %default]')
    parser.add_option('--orthologs', dest='orthologMode', default=False, action='store_true', help='Mouse file to prepare human and mouse 1-1 ortholog set (Mouse_FantomAnnotations.InputData.pM10Kb.txt.gz recommended to accompany human data_file): [Default: %default]')
    parser.add_option('--over', dest='overwrite', default=False, action='store_true', help='Overwrite directory? [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        print(args)
        parser.error('Must provide data file and output directory')
    else:
        data_file = args[0]
        out_dir = args[1]
        compress_args = {'compression': 'gzip', 'compression_opts': 1}
        trainfile = os.path.join(out_dir, 'train.h5')
        validfile = os.path.join(out_dir, 'valid.h5')
        testfile = os.path.join(out_dir, 'test.h5')

    if options.orthologMode:
        trainfile = os.path.join(out_dir, 'train_human1to1.h5')
        validfile = os.path.join(out_dir, 'valid_human1to1.h5')
        testfile = os.path.join(out_dir, 'test_human1to1.h5')
        trainfile2 = os.path.join(out_dir, 'train_mouse1to1.h5')
        validfile2 = os.path.join(out_dir, 'valid_mouse1to1.h5')
        testfile2 = os.path.join(out_dir, 'test_mouse1to1.h5')

    if options.overwrite or not os.path.exists(out_dir):
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        # load data
        promoters, halflifedata, labels, geneNames = preprocess(data_file, options.orthologMode)

        # check that the sum is valid
        assert(options.testCount + options.validCount <= promoters.shape[0])
        test_count = options.testCount
        valid_count = options.validCount

        train_count = promoters.shape[0] - test_count - valid_count

        if options.crossVal:
            print('running 10-fold cross val w/ %d sequences ' % promoters.shape[0])
            kf = KFold(n_splits=10, random_state=42, shuffle=False)
            fold = 0
            for train_index, test_index in kf.split(promoters): #keep aside 1000 examples of train indices for validation set
                fold += 1
                print('fold %d' % fold)
                h5f_train = h5py.File(os.path.join(out_dir, str(fold)+'train.h5'), 'w')
                h5f_valid = h5py.File(os.path.join(out_dir, str(fold)+'valid.h5'), 'w')
                h5f_test = h5py.File(os.path.join(out_dir, str(fold)+'test.h5'), 'w')
                valid_index = train_index[0:1000]
                train_index = train_index[1000:len(train_index)]
                h5f_train.create_dataset('data'    , data=halflifedata[train_index,:], **compress_args)
                h5f_train.create_dataset('promoter', data=promoters[train_index,:], **compress_args)
                h5f_train.create_dataset('label'   , data=labels[train_index], **compress_args)
                h5f_train.create_dataset('geneName' , data=np.array(geneNames)[train_index].tolist(), **compress_args)
                h5f_train.close()
                h5f_valid.create_dataset('data'    , data=halflifedata[valid_index,:], **compress_args)
                h5f_valid.create_dataset('promoter', data=promoters[valid_index,:], **compress_args)
                h5f_valid.create_dataset('label'   , data=labels[valid_index], **compress_args)
                h5f_valid.create_dataset('geneName' , data=np.array(geneNames)[valid_index].tolist(), **compress_args)
                h5f_valid.close()
                h5f_test.create_dataset('data'    , data=halflifedata[test_index,:], **compress_args)
                h5f_test.create_dataset('promoter', data=promoters[test_index,:], **compress_args)
                h5f_test.create_dataset('label'   , data=labels[test_index], **compress_args)
                h5f_test.create_dataset('geneName' , data=np.array(geneNames)[test_index].tolist(), **compress_args)
                h5f_test.close()
        else:
            print('%d training sequences ' % train_count)
            print('%d test sequences ' % test_count)
            print('%d validation sequences ' % valid_count)
            h5f_train = h5py.File(trainfile, 'w')
            h5f_valid = h5py.File(validfile, 'w')
            h5f_test = h5py.File(testfile, 'w')
            i = 0
            if train_count > 0:
                h5f_train.create_dataset('data'    , data=halflifedata[i:i+train_count,:], **compress_args)
                h5f_train.create_dataset('promoter', data=promoters[i:i+train_count,:], **compress_args)
                h5f_train.create_dataset('label'   , data=labels[i:i+train_count], **compress_args)
                h5f_train.create_dataset('geneName' , data=geneNames[i:i+train_count], **compress_args)
                h5f_train.close()
            i += train_count
            if valid_count > 0:
                h5f_valid.create_dataset('data'    , data=halflifedata[i:i+valid_count,:], **compress_args)
                h5f_valid.create_dataset('promoter', data=promoters[i:i+valid_count,:], **compress_args)
                h5f_valid.create_dataset('label'   , data=labels[i:i+valid_count], **compress_args)
                h5f_valid.create_dataset('geneName' , data=geneNames[i:i+valid_count], **compress_args)
                h5f_valid.close()
            i += valid_count
            if test_count > 0:
                h5f_test.create_dataset('data'    , data=halflifedata[i:i+test_count,:], **compress_args)
                h5f_test.create_dataset('promoter', data=promoters[i:i+test_count,:], **compress_args)
                h5f_test.create_dataset('label'   , data=labels[i:i+test_count], **compress_args)
                h5f_test.create_dataset('geneName' , data=geneNames[i:i+test_count], **compress_args)
                h5f_test.close()

            if options.orthologMode:
                print("Finding 1-1 orthologs...")
                h5f_train = h5py.File(trainfile2, 'w')
                h5f_valid = h5py.File(validfile2, 'w')
                h5f_test = h5py.File(testfile2, 'w')
                promoters2, halflifedata2, labels2, geneNames2 = preprocess(options.orthologMode, options.orthologMode)
                orthologs = pd.read_table("human2mouse_one2one_orthologs.txt", header=None)

                i = 0
                orthoids = orthologs[orthologs[0].isin(geneNames[i:i+train_count])][1] #transform human to mouse IDs
                idxs = np.where(np.isin(geneNames2,orthoids))[0].tolist()
                h5f_train.create_dataset('data'    , data=halflifedata2[idxs,:], **compress_args)
                h5f_train.create_dataset('promoter', data=promoters2[idxs,:], **compress_args)
                h5f_train.create_dataset('label'   , data=labels2[idxs], **compress_args)
                h5f_train.create_dataset('geneName', data=np.array(geneNames2)[idxs].tolist(), **compress_args)
                print('%d 1-1 mouse orthologs found for training set' % labels2[idxs].shape)
                h5f_train.close()
                i += train_count
                orthoids = orthologs[orthologs[0].isin(geneNames[i:i+valid_count])][1]
                idxs = np.where(np.isin(geneNames2,orthoids))[0].tolist()
                h5f_valid.create_dataset('data'    , data=halflifedata2[idxs,:], **compress_args)
                h5f_valid.create_dataset('promoter', data=promoters2[idxs,:], **compress_args)
                h5f_valid.create_dataset('label'   , data=labels2[idxs], **compress_args)
                h5f_valid.create_dataset('geneName', data=np.array(geneNames2)[idxs].tolist(), **compress_args)
                print('%d 1-1 mouse orthologs found validation set' % labels2[idxs].shape)
                h5f_valid.close()
                i += valid_count
                orthoids = orthologs[orthologs[0].isin(geneNames[i:i+test_count])][1]
                idxs = np.isin(geneNames2,orthoids)
                h5f_test.create_dataset('data'    , data=halflifedata2[idxs,:], **compress_args)
                h5f_test.create_dataset('promoter', data=promoters2[idxs,:], **compress_args)
                h5f_test.create_dataset('label'   , data=labels2[idxs], **compress_args)
                h5f_test.create_dataset('geneName', data=np.array(geneNames2)[idxs].tolist(), **compress_args)
                print('%d 1-1 mouse orthologs found for test set ' % labels2[idxs].shape)
                h5f_test.close()
    else:
        parser.error('Nothing done...Run with --over to overwrite')

def one_hot(seq):
    seq_len = len(seq.item(0))
    seqindex = {'A':0, 'C':1, 'G':2, 'T':3, 'a':0, 'c':1, 'g':2, 't':3}
    seq_vec = np.zeros((len(seq),seq_len,4), dtype='bool')
    for i in range(len(seq)):
        thisseq = seq.item(i)
        for j in range(seq_len):
            try:
                seq_vec[i,j,seqindex[thisseq[j]]] = 1
            except:
                pass
    return seq_vec

def preprocess(data_file, orthologMode):
    table = pd.read_table(data_file, index_col=0)
    maskedIDs = pd.read_table("mask_histone_genes_mm10.txt", header=None) #mask histone genes, chrY genes already filtered out
    maskedIDs2 = pd.read_table("mask_histone_genes.txt", header=None) #mask histone genes, chrY genes already filtered out
    table = table[~table.index.isin(maskedIDs[0])] #remove rows corresponding to chrY or histone sequences
    table = table[~table.index.isin(maskedIDs2[0])] #remove rows corresponding to chrY or histone sequences
    if orthologMode:
        orthologs = pd.read_table("1to1_orthologs_expression.txt", header=None)
        table = table[table.index.isin(orthologs[[0,1]].values.flatten())] #must match human or mouse 1-1 ortholog IDs
    table[table.columns[range(0,5)+[8]]] = np.log10(table[table.columns[range(0,5)+[8]]]+0.1)
    table = table.sample(table.shape[0], replace=False, random_state=1)
    table[table.columns[range(0,9)]] = preprocessing.scale(table[table.columns[range(0,9)]])
    print("\nPre-processed data...one-hot encoding...")
    promoters = one_hot(table['PROMOTER'].as_matrix())
    halflifedata = table[table.columns[range(1,9)]].as_matrix()
    labels = table['EXPRESSION'].as_matrix()
    geneNames = list(table.index)
    print("Processed data from %s" % data_file)
    return promoters, halflifedata, labels, geneNames

if __name__ == '__main__':
    main()
