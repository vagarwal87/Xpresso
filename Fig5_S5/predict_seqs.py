import sys, pickle
import pandas as pd
import numpy as np
from optparse import OptionParser
from keras.models import Model, load_model

def main():
    usage = 'usage: %prog [options] <param_file> <trained_model> <test_file> <outfile>'
    parser = OptionParser(usage)
    parser.add_option('--revCom', dest='revcom', default=False, action='store_true', help='Make predictions for minus strand instead of plus? % [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 4:
        print args
        parser.error('Must provide mode hyperparameter file and 2-column file to generate predictions for')
    else:
        param_file = args[0]
        trained_model = args[1]
        test_file = args[2]
        outfile = args[3]

    def revCom(x):
        for y in range(0,x.shape[0]):
            x[y] = np.fliplr(np.flipud(x[y]))
        return x

    trials = pickle.load(open(param_file, "rb"))
    best = trials.argmin
    model = load_model(trained_model, compile=False)

    table = pd.read_table(test_file, index_col=0, header=None)
    seqs = one_hot(table.as_matrix())
    if options.revcom:
        seqs = revCom(seqs)
    if seqs.shape[1] != 10500:
        tsspos = 7000
        leftpos = tsspos - seqs.shape[1] / 2
        if seqs.shape[1] <= tsspos:
            tmpseqs = np.zeros((seqs.shape[0],10500,4), dtype='bool')
            tmpseqs[:,leftpos:(leftpos+seqs.shape[1]),:] = seqs
            seqs = tmpseqs
        else:
            print >> sys.stderr, 'Sequences are above the allowable size of 10500nt'
            sys.exit()
    halflifedata = np.zeros((seqs.shape[0],6), dtype='float16')
    print >> sys.stderr, "Processed data from %s" % test_file
    predictions_test = model.predict([seqs, halflifedata], batch_size=20).flatten()
    df = pd.DataFrame(np.column_stack((table.index, predictions_test)), columns=['Info','Pred'])
    df.to_csv(outfile, index=False, header=True, sep='\t')

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

if __name__ == '__main__':
    main()
