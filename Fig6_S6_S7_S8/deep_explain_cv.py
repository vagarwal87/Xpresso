import os, h5py, math
import numpy as np
import pandas as pd
import tensorflow as tf
from optparse import OptionParser
from tensorflow import keras
from keras.models import load_model, Model
from keras import backend as K
from deepexplain.tensorflow import DeepExplain

batchsize = 500

def main():
    usage = 'usage: %prog [options] <data_file> <out_dir> <cv_fold>'
    parser = OptionParser(usage)
    (options,args) = parser.parse_args()

    if len(args) != 3:
        print args
        parser.error('Must provide data file and output directory')
    else:
        data_file = args[0]
        out_dir = args[1]
        fold = args[2]
        testfile = os.path.join(out_dir, fold+'test.h5')

    testfile = h5py.File(testfile, 'r')
    X_testhalflife, X_testpromoter, y_test, geneName = testfile['data'], testfile['promoter'], testfile['label'], testfile['geneName']
    model = load_model(data_file)

    with DeepExplain(session=K.get_session()) as de:
        input_tensor = model.inputs
        fModel = Model(inputs = input_tensor, outputs = model.outputs)
        for method in ['deeplift', 'grad*input', 'saliency', 'elrp', 'intgrad']: #'occlusion' not supported
            pdframe = pd.DataFrame()
            for i in range(0, int(math.ceil(len(geneName) / float(batchsize)))):
                first = i*batchsize
                last = (i*batchsize+batchsize)
                if last > len(geneName): last = len(geneName)
                xs = X_testpromoter[first:last,3000:13500,:]
                xs2 = X_testhalflife[first:last,:]
                ys = y_test[first:last]
                gN = geneName[first:last]
                if method in ('intgrad', 'deeplift'): #try these methods with and without specified baseline
                    #empirical ACGT frequencies in -7Kb to +3.5Kb sequence surrounding human TSSs (for non-expressed genes only)
                    baseline = [np.repeat(np.array([[0.2617064, 0.2335449, 0.2379253, 0.2668234]]), 10500, axis=0), np.zeros(6)]
                    map = de.explain(method, fModel(input_tensor), input_tensor, [xs, xs2], baseline = baseline)
                else:
                    map = de.explain(method, fModel(input_tensor), input_tensor, [xs, xs2])
                X, halflife = map[0], map[1]
                frame = pd.DataFrame(np.column_stack((ys, 10**3 * np.sum(X, 2), halflife)))
                frame.index = gN
                pdframe = pdframe.append(frame)
            pdframe.to_csv(os.path.join(out_dir, method.replace("*", "")+'.'+fold+'.txt'),sep='\t',header=False, float_format='%.3f')

if __name__ == '__main__':
    main()
