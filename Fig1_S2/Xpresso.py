import sys, os, h5py, pickle
import pandas as pd
from optparse import OptionParser
from scipy import stats
import tensorflow as tf
from tensorflow import keras
from keras.optimizers import Adam
from keras.models import Model, load_model
from keras.layers import *
from keras.metrics import *
from keras.utils import plot_model
from keras import backend as K
from keras.callbacks import Callback, ModelCheckpoint, EarlyStopping
from hyperopt import fmin, tpe, rand, anneal, hp, STATUS_OK, STATUS_FAIL, Trials, mix, partial, space_eval

global X_trainhalflife, X_trainpromoter, y_train, geneName_train, X_validhalflife, X_validpromoter, y_valid, geneName_valid, X_testhalflife, X_testpromoter, y_test, geneName_test, params

def main():
    usage = 'usage: %prog [options] <mode> <database_file> <database_directory>'
    parser = OptionParser(usage)
    parser.add_option('-c', dest='counts', default=0, type='int', help='Number of training counts to subsample [Default: %default]')
    parser.add_option('--bestmanual', dest='bestmanual', default=False, action='store_true', help='Try best manually identified model % [Default: %default]')
    parser.add_option('--fold', dest='cvfold', default='', type='string', help='Which of the 10 folds of cross-validation to use % [Default: %default]')
    parser.add_option('--trial', dest='trial', default='', type='string', help='Trial number % [Default: %default]')
    parser.add_option('--usemodel', dest='usemodel', default=None, type='string', help='Use pre-trained model % [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 3:
        print (args)
        parser.error('Must provide mode (tune, train, or test), hyperparameter database file, and database directory')
    else:
        mode = args[0]
        database = args[1]
        datadir = args[2]

    global X_trainhalflife, X_trainpromoter, y_train, geneName_train, X_validhalflife, X_validpromoter, y_valid, geneName_valid, X_testhalflife, X_testpromoter, y_test, geneName_test, params
    params['datadir'] = datadir
    if not options.usemodel:
        trainfile = h5py.File(os.path.join(datadir, options.cvfold+'train.h5'), 'r') #_mouse1to1
        X_trainhalflife, X_trainpromoter, y_train, geneName_train = trainfile['data'], trainfile['promoter'], trainfile['label'], trainfile['geneName']
        validfile = h5py.File(os.path.join(datadir, options.cvfold+'valid.h5'), 'r') #_mouse1to1
        X_validhalflife, X_validpromoter, y_valid, geneName_valid = validfile['data'], validfile['promoter'], validfile['label'], validfile['geneName']

    if mode == "tune":
        while True: # loop indefinitely and stop whenever you like
            run_trials(database)
    else:
        testfile = h5py.File(os.path.join(datadir, options.cvfold+'test.h5'), 'r') #_mouse1to1_human1to1
        X_testhalflife, X_testpromoter, y_test, geneName_test = testfile['data'], testfile['promoter'], testfile['label'], testfile['geneName']
        if options.bestmanual:
            params = { 'datadir' : datadir, 'batchsize' : 2**6, 'leftpos' : 8500, 'rightpos' : 11500, 'activationFxn' : 'relu', 'numFiltersConv1' : 2**6, 'filterLenConv1' : 5, 'dilRate1' : 1,
                       'maxPool1' : 10, 'numconvlayers' : { 'numFiltersConv2' : 2**6, 'filterLenConv2' : 5, 'dilRate2' : 1, 'maxPool2' : 20, 'numconvlayers1' : { 'numconvlayers2' : 'two' } },
                       'dense1' : 100, 'dropout1' : 0.5, 'numdenselayers' : { 'layers' : 'one' } }
            print("Using best human-identified parameters")
        else:
            trials = pickle.load(open(database, "rb"))
            best = trials.argmin
            params = space_eval(params, best)
            print("Found saved Trials!")
        print ("The best parameters are:")
        print (params)
        params['subsample'] = options.counts
        params['cvfold'] = options.cvfold
        params['trial'] = options.trial
        params['usemodel'] = options.usemodel
        params['tuneMode'] = 0 #enable mode that trains best model structure over up to 100 epochs, and evaluates final model on test set
        results = objective(params)
        print("Best Validation MSE = %.3f" % results['loss'])

params = {
    'tuneMode' : 1,
    'batchsize' : 2**hp.quniform('batchsize', 5, 7, 1),
    'leftpos' : hp.quniform('leftpos', 0, 10000, 500),
    'rightpos' : hp.quniform('rightpos', 10000, 20000, 500),
    'activationFxn' : 'relu', #hp.choice('activationFxn', ['relu', 'elu', 'selu', 'LeakyReLU', 'PReLU']) -- tried but none worked better than simply relu
    'numFiltersConv1' : 2**hp.quniform('numFiltersConv1', 4, 7, 1),
    'filterLenConv1' : hp.quniform('filterLenConv1', 1, 10, 1),
    'dilRate1' : hp.quniform('dilRate1', 1, 4, 1),
    'maxPool1' : hp.quniform('maxPool1', 5, 100, 5),
    'numconvlayers' : hp.choice('numconvlayers', [
    {
        'numconvlayers1' : 'one'
    },
    {
        'numFiltersConv2' : 2**hp.quniform('numFiltersConv2', 4, 7, 1),
        'filterLenConv2' : hp.quniform('filterLenConv2', 1, 10, 1),
        'dilRate2' : hp.quniform('dilRate2', 1, 4, 1),
        'maxPool2' : hp.quniform('maxPool2', 5, 100, 5),
        'numconvlayers1' : hp.choice('numconvlayers1', [
        {
            'numconvlayers2' : 'two'
        },
        {
            'numFiltersConv3' : 2**hp.quniform('numFiltersConv3', 4, 7, 1),
            'filterLenConv3' : hp.quniform('filterLenConv3', 1, 10, 1),
            'dilRate3' : hp.quniform('dilRate3', 1, 4, 1),
            'maxPool3' : hp.quniform('maxPool3', 5, 100, 5),
            'numconvlayers2' : hp.choice('numconvlayers2', [
            {
                'numconvlayers3' : 'three'
            },
            {
                'numFiltersConv4' : 2**hp.quniform('numFiltersConv4', 4, 7, 1),
                'filterLenConv4' : hp.quniform('filterLenConv4', 1, 10, 1),
                'dilRate4' : hp.quniform('dilRate4', 1, 4, 1),
                'maxPool4' : hp.quniform('maxPool4', 5, 100, 5),
                'numconvlayers3' : 'four'
            }])
        }])
    }]),
    'dense1' : 2**hp.quniform('dense1', 1, 8, 1),
    'dropout1' : hp.uniform('dropout1', 0, 1),
    'numdenselayers' : hp.choice('numdenselayers', [
        {
            'layers' : 'one'
        },
        {
            'layers' : 'two' ,
            'dense2' : 2**hp.quniform('dense2', 1, 8, 1),
            'dropout2' : hp.uniform('dropout2', 0, 1)
        }
    ])
}

def run_trials(database):
    trials_step = 5  # how many additional trials to do after loading saved trials
    max_trials = 5  # initial max_trials. put something small to not have to wait

    try:  # try to load an already saved trials object, and increase the max
        trials = pickle.load(open(database, "rb"))
        print("Found saved Trials! Loading...")
        max_trials = len(trials.trials) + trials_step
        print("Rerunning from {} trials to {} (+{}) trials".format(len(trials.trials), max_trials, trials_step))
    except:  # create a new trials object and start searching
        trials = Trials()

    best = fmin(objective, params, max_evals = max_trials, trials = trials,
        algo = anneal.suggest)
        # algo = rand.suggest)
        # algo = tpe.suggest)
        # algo = partial(mix.suggest, p_suggest=[(0.2, rand.suggest),(0.6, tpe.suggest),(0.2, anneal.suggest)]))

    ##### sample random parameter sets and print
    # import hyperopt.pyll.stochastic
    # print (hyperopt.pyll.stochastic.sample(params))

    print( "Best:", best)
    # save the trials object
    with open(database, "wb") as f:
        pickle.dump(trials, f)

def objective(params):
    leftpos = int(params['leftpos'])
    rightpos = int(params['rightpos'])
    activationFxn = params['activationFxn']
    if not params['usemodel']:
        global X_trainhalflife, y_train
        X_trainpromoterSubseq = X_trainpromoter[:,leftpos:rightpos,:]
        X_validpromoterSubseq = X_validpromoter[:,leftpos:rightpos,:]
        halflifedata = Input(shape=(X_trainhalflife.shape[1:]), name='halflife')
        input_promoter = Input(shape=X_trainpromoterSubseq.shape[1:], name='promoter')

    try:
    # if True:
        mse = 1
        if params['usemodel']:
            model = load_model(params['usemodel'])
            print('Loaded results from:', params['usemodel'])
        else:
            x = Conv1D(int(params['numFiltersConv1']), int(params['filterLenConv1']), dilation_rate=int(params['dilRate1']), padding='same', kernel_initializer='glorot_normal', input_shape=X_trainpromoterSubseq.shape[1:],activation=activationFxn)(input_promoter)
            x = MaxPooling1D(int(params['maxPool1']))(x)

            if params['numconvlayers']['numconvlayers1'] != 'one':
                maxPool2 = int(params['numconvlayers']['maxPool2'])
                x = Conv1D(int(params['numconvlayers']['numFiltersConv2']), int(params['numconvlayers']['filterLenConv2']), dilation_rate=int(params['numconvlayers']['dilRate2']), padding='same', kernel_initializer='glorot_normal',activation=activationFxn)(x) #[2, 3, 4, 5, 6, 7, 8, 9, 10]
                x = MaxPooling1D(maxPool2)(x)
                if params['numconvlayers']['numconvlayers1']['numconvlayers2'] != 'two':
                    maxPool3 = int(params['numconvlayers']['numconvlayers1']['maxPool3'])
                    x = Conv1D(int(params['numconvlayers']['numconvlayers1']['numFiltersConv3']), int(params['numconvlayers']['numconvlayers1']['filterLenConv3']), dilation_rate=int(params['numconvlayers']['numconvlayers1']['dilRate3']), padding='same', kernel_initializer='glorot_normal',activation=activationFxn)(x) #[2, 3, 4, 5]
                    x = MaxPooling1D(maxPool3)(x)
                    if params['numconvlayers']['numconvlayers1']['numconvlayers2']['numconvlayers3'] != 'three':
                        maxPool4 = int(params['numconvlayers']['numconvlayers1']['numconvlayers2']['maxPool4'])
                        x = Conv1D(int(params['numconvlayers']['numconvlayers1']['numconvlayers2']['numFiltersConv4']), int(params['numconvlayers']['numconvlayers1']['numconvlayers2']['filterLenConv4']), dilation_rate=int(params['numconvlayers']['numconvlayers1']['numconvlayers2']['dilRate4']), padding='same', kernel_initializer='glorot_normal',activation=activationFxn)(x) #[2, 3, 4, 5]
                        x = MaxPooling1D(maxPool4)(x)

            x = Flatten()(x)
            x = Concatenate()([x, halflifedata])
            x = Dense(int(params['dense1']))(x)
            x = Activation(activationFxn)(x)
            x = Dropout(params['dropout1'])(x)
            if params['numdenselayers']['layers'] == 'two':
                x = Dense(int(params['numdenselayers']['dense2']))(x)
                x = Activation(activationFxn)(x)
                x = Dropout(params['numdenselayers']['dropout2'])(x)
            main_output = Dense(1)(x)
            model = Model(inputs=[input_promoter, halflifedata], outputs=[main_output])
            model.compile(Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay=0.0),'mean_squared_error', metrics=['mean_squared_error'])

        if params['tuneMode']:
            result = model.fit([X_trainpromoterSubseq, X_trainhalflife], y_train, batch_size=int(params['batchsize']), shuffle="batch", epochs=10,
                                validation_data=[[X_validpromoterSubseq, X_validhalflife], y_valid])
            mse = min(result.history['val_mean_squared_error'])
            print("leftpos, rightpos, mse")
            print(leftpos, rightpos, mse)
        else:
            print(model.summary())
            plot_model(model, to_file=os.path.join(params['datadir'], 'best_model.png')) #requires Pydot/Graphviz to generate graph of network
            X_testpromoterSubseq = X_testpromoter[:,leftpos:rightpos,:]
            if not params['usemodel']:
                if params['subsample'] > 0:
                    X_trainpromoterSubseq = X_trainpromoterSubseq[0:params['subsample'],:,:]
                    X_trainhalflife = X_trainhalflife[0:params['subsample'],:]
                    y_train = y_train[0:params['subsample']]
                check_cb = ModelCheckpoint(os.path.join(params['datadir'], params['trial']+params['cvfold']+'trainepoch.{epoch:02d}-{val_loss:.4f}.h5'), monitor='val_loss', verbose=1, save_best_only=True, mode='min')
                earlystop_cb = EarlyStopping(monitor='val_loss', patience=5, verbose=1, mode='min')
                result = model.fit([X_trainpromoterSubseq, X_trainhalflife], y_train, batch_size=int(params['batchsize']), shuffle="batch", epochs=100,
                    validation_data=[[X_validpromoterSubseq, X_validhalflife], y_valid], callbacks=[earlystop_cb, check_cb])
                mse_history = result.history['val_mean_squared_error']
                mse = min(mse_history)
                best_file = os.path.join(params['datadir'], params['trial']+params['cvfold']+'trainepoch.%02d-%.4f.h5' % (mse_history.index(mse), mse))
                model = load_model(best_file)
                print('Loaded results from:', best_file)

            predictions_test = model.predict([X_testpromoterSubseq, X_testhalflife], batch_size=20).flatten()
            slope, intercept, r_value, p_value, std_err = stats.linregress(predictions_test, y_test)
            print('Test R^2 = %.3f' % r_value**2)
            df = pd.DataFrame(np.column_stack((geneName_test, predictions_test, y_test)), columns=['Gene','Pred','Actual'])
            df.to_csv(os.path.join(params['datadir'], params['trial']+params['cvfold']+'predictions.txt'), index=False, header=True, sep='\t')

        return {'loss': mse, 'status': STATUS_OK }

    except:
        return {'loss': 1, 'status': STATUS_FAIL } # loss = 1 indicates a poor-performing model; reason model might fail include: incompatible parameters or insufficient memory resources available

if __name__ == '__main__':
    main()
