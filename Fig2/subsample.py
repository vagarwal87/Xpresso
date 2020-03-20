import sys, os, h5py
from numpy.random import choice
import numpy as np
from optparse import OptionParser

def main():
    usage = 'usage: %prog [options] <subsamples> <indir> <outdir>'
    parser = OptionParser(usage)
    (options,args) = parser.parse_args()
    if len(args) != 3:
        print(args)
        sys.exit('Must provide number of subsamples, input directory, and output directory')
    counts, in_dir, out_dir = args
    counts = int(counts)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    file = h5py.File(os.path.join(in_dir, 'train.h5'), 'r')
    X_halflife, X_promoter, y, geneName = file['data'], file['promoter'], file['label'], file['geneName']
    for i in range(10):
        print('sample %d' % i)
        os.symlink(os.path.join(os.path.realpath(in_dir), 'valid.h5'), os.path.join(out_dir, str(i)+'valid.h5'))
        os.symlink(os.path.join(os.path.realpath(in_dir), 'test.h5'), os.path.join(out_dir, str(i)+'test.h5'))
        h5f = h5py.File(os.path.join(out_dir, str(i)+'train.h5'), 'w')
        index = np.sort(choice(X_promoter.shape[0], size=counts, replace=False))
        compress_args = {'compression': 'gzip', 'compression_opts': 1}
        h5f.create_dataset('data'    , data=X_halflife[index,:], **compress_args)
        h5f.create_dataset('promoter', data=X_promoter[index,:,:], **compress_args)
        h5f.create_dataset('label'   , data=y[index.tolist()], **compress_args)
        h5f.create_dataset('geneName', data=geneName[index.tolist()], **compress_args)
        h5f.close()

if __name__ == '__main__':
    main()
