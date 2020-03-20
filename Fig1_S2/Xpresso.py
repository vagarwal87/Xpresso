import sys, os, h5py, pickle
from optparse import OptionParser
from scipy import stats
from keras.optimizers import Adam
from keras.models import Model, load_model
from keras.layers import *
from keras.metrics import *
from keras import backend as K
from keras.callbacks import Callback, ModelCheckpoint, EarlyStopping
from hyperopt import fmin, tpe, rand, anneal, hp, STATUS_OK, STATUS_FAIL, Trials, mix, partial, space_eval

global X_trainhalflife, X_trainpromoter, y_train, headers_train, X_validhalflife, X_validpromoter, y_valid, headers_valid, X_testhalflife, X_testpromoter, y_test, headers_test, params

def main():
    usage = 'usage: %prog [options] <train> <database_file> <database_directory>'
    parser = OptionParser(usage)
    (options,args) = parser.parse_args()

    if len(args) != 3:
        print args
        parser.error('Must provide mode (train or test), hyperparameter database file, and database directory')
    else:
        train = args[0]
        database = args[1]
        datadir = args[2]

    global X_trainhalflife, X_trainpromoter, y_train, headers_train, X_validhalflife, X_validpromoter, y_valid, headers_valid, X_testhalflife, X_testpromoter, y_test, headers_test, params
    params['datadir'] = datadir
    trainfile = h5py.File(os.path.join(datadir, 'train.h5'), 'r')
    X_trainhalflife, X_trainpromoter, y_train, headers_train = trainfile['data'], trainfile['promoter'], trainfile['label'], trainfile['headers']
    validfile = h5py.File(os.path.join(datadir, 'valid.h5'), 'r')
    X_validhalflife, X_validpromoter, y_valid, headers_valid = validfile['data'], validfile['promoter'], validfile['label'], validfile['headers']

    if train == "train":
        while True: # loop indefinitely and stop whenever you like
            run_trials(database)
    elif train == "test":
        testfile = h5py.File(os.path.join(datadir, 'test.h5'), 'r')
        X_testhalflife, X_testpromoter, y_test, headers_test = testfile['data'], testfile['promoter'], testfile['label'], testfile['headers']
        params['trainMode'] = 0 #enable mode that trains best model structure over up to 100 epochs, and evaluates final model on test set
        # trials = pickle.load(open(database, "rb"))
        # best = trials.argmin
        print "Found saved Trials! The best parameters are:"
        # params = space_eval(params, best)
        # print params
        human_best = { 'datadir' : datadir, 'trainMode' : 0, 'batchsize' : 2**5, 'leftpos' : 8500, 'rightpos' : 11500, 'activationFxn' : 'relu', 'numFiltersConv1' : 2**6, 'filterLenConv1' : 5, 'dilRate1' : 1,
                   'maxPool1' : 20, 'numconvlayers' : { 'numFiltersConv2' : 2**6, 'filterLenConv2' : 5, 'dilRate2' : 1, 'maxPool2' : 10, 'numconvlayers1' : { 'numconvlayers2' : 'two' } },
                   'dense1' : 100, 'dropout1' : 0.5, 'numdenselayers' : { 'layers' : 'one' } }
        results = objective(human_best)
        print "Best Validation MSE = %.3f", results['loss']

params = {
    'trainMode' : 1,
    'batchsize' : 2**hp.quniform('batchsize', 5, 7, 1),
    'leftpos' : hp.quniform('leftpos', 0, 10000, 200),
    'rightpos' : hp.quniform('rightpos', 10000, 20000, 200),
    'activationFxn' : hp.choice('activationFxn', ['relu', 'elu', 'selu', 'tanh']),
    'numFiltersConv1' : 2**hp.quniform('numFiltersConv1', 4, 7, 1),
    'filterLenConv1' : hp.quniform('filterLenConv1', 1, 20, 1),
    'dilRate1' : hp.quniform('dilRate1', 1, 5, 1),
    'maxPool1' : hp.quniform('maxPool1', 1, 20, 1),
    'numconvlayers' : hp.choice('numconvlayers', [
    {
        'numconvlayers1' : 'one'
    },
    {
        'numFiltersConv2' : 2**hp.quniform('numFiltersConv2', 4, 7, 1),
        'filterLenConv2' : hp.quniform('filterLenConv2', 1, 20, 1),
        'dilRate2' : hp.quniform('dilRate2', 1, 5, 1),
        'maxPool2' : hp.quniform('maxPool2', 1, 20, 1),
        'numconvlayers1' : hp.choice('numconvlayers1', [
        {
            'numconvlayers2' : 'two'
        },
        {
            'numFiltersConv3' : 2**hp.quniform('numFiltersConv3', 4, 7, 1),
            'filterLenConv3' : hp.quniform('filterLenConv3', 1, 20, 1),
            'dilRate3' : hp.quniform('dilRate3', 1, 5, 1),
            'maxPool3' : hp.quniform('maxPool3', 1, 20, 1),
            'numconvlayers2' : hp.choice('numconvlayers2', [
            {
                'numconvlayers3' : 'three'
            },
            {
                'numFiltersConv4' : 2**hp.quniform('numFiltersConv4', 4, 7, 1),
                'filterLenConv4' : hp.quniform('filterLenConv4', 1, 20, 1),
                'dilRate4' : hp.quniform('dilRate4', 1, 5, 1),
                'maxPool4' : hp.quniform('maxPool4', 1, 20, 1),
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
    trials_step = 5  # how many additional trials to do after loading saved trials. 1 = save after iteration
    max_trials = 5  # initial max_trials. put something small to not have to wait

    try:  # try to load an already saved trials object, and increase the max
        trials = pickle.load(open(database, "rb"))
        print "Found saved Trials! Loading..."
        max_trials = len(trials.trials) + trials_step
        print "Rerunning from {} trials to {} (+{}) trials".format(len(trials.trials), max_trials, trials_step)
    except:  # create a new trials object and start searching
        trials = Trials()

    best = fmin(objective, params, max_evals = max_trials, trials = trials,
        algo = partial(mix.suggest, p_suggest=[(0.2, rand.suggest),(0.6, tpe.suggest),(0.2, anneal.suggest)])) #tpe.suggest

    ##### sample random parameter sets and print
    # import hyperopt.pyll.stochastic
    # print hyperopt.pyll.stochastic.sample(params)

    print "Best:", best
    # save the trials object
    with open(database, "wb") as f:
        pickle.dump(trials, f)

def objective(params):
    leftpos = int(params['leftpos'])
    rightpos = int(params['rightpos'])
    activationFxn = params['activationFxn']

    X_trainpromoterSubseq = X_trainpromoter[:,leftpos:rightpos,:]
    X_validpromoterSubseq = X_validpromoter[:,leftpos:rightpos,:]

    halflifedata = Input(shape=(X_trainhalflife.shape[1:]), name='halflife')
    input_promoter = Input(shape=X_trainpromoterSubseq.shape[1:], name='promoter')
    # runFail =

    # try:
    x = Conv1D(int(params['numFiltersConv1']), int(params['filterLenConv1']), dilation_rate=int(params['dilRate1']), padding='same', kernel_initializer='glorot_normal', input_shape=X_trainpromoterSubseq.shape[1:],activation=activationFxn)(input_promoter)
    # BatchNormalization()(x)
    maxPool1 = int(params['maxPool1'])
    x = MaxPooling1D()(x)

    if params['numconvlayers']['numconvlayers1'] != 'one':
        maxPool2 = int(params['numconvlayers']['maxPool2'])
        # if (maxPool1 * maxPool2 > rightpos - leftpos): return runFail
        x = Conv1D(int(params['numconvlayers']['numFiltersConv2']), int(params['numconvlayers']['filterLenConv2']), dilation_rate=int(params['numconvlayers']['dilRate2']), padding='same', kernel_initializer='glorot_normal',activation=activationFxn)(x) #[2, 3, 4, 5, 6, 7, 8, 9, 10]
        x = MaxPooling1D(maxPool2)(x)
        if params['numconvlayers']['numconvlayers1']['numconvlayers2'] != 'two':
            maxPool3 = int(params['numconvlayers']['numconvlayers1']['maxPool3'])
            # if (maxPool1 * maxPool2 * maxPool3 > rightpos - leftpos): return runFail
            x = Conv1D(int(params['numconvlayers']['numconvlayers1']['numFiltersConv3']), int(params['numconvlayers']['numconvlayers1']['filterLenConv3']), dilation_rate=int(params['numconvlayers']['numconvlayers1']['dilRate3']), padding='same', kernel_initializer='glorot_normal',activation=activationFxn)(x) #[2, 3, 4, 5]
            x = MaxPooling1D(maxPool3)(x)
            if params['numconvlayers']['numconvlayers1']['numconvlayers2']['numconvlayers3'] != 'three':
                maxPool4 = int(params['numconvlayers']['numconvlayers1']['numconvlayers2']['maxPool4'])
                # if (maxPool1 * maxPool2 * maxPool3 * maxPool4 > rightpos - leftpos): return runFail
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

    if params['trainMode']:
        result = model.fit([X_trainpromoterSubseq, X_trainhalflife], y_train, batch_size=int(params['batchsize']), shuffle="batch", epochs=1,
                            validation_data=[[X_validpromoterSubseq, X_validhalflife], y_valid])
        mse = result.history['val_mean_squared_error'][-1]
        print "leftpos, rightpos, mse"
        print leftpos, rightpos, mse
    else:
        X_testpromoterSubseq = X_testpromoter[:,leftpos:rightpos,:]
        check_cb = ModelCheckpoint(os.path.join(params['datadir'], 'trainepoch.{epoch:02d}.h5'), monitor='val_loss', verbose=1, save_best_only=True, mode='min')
        earlystop_cb = EarlyStopping(monitor='val_loss', patience=7, verbose=1, mode='min')
        result = model.fit([X_trainpromoterSubseq, X_trainhalflife], y_train, batch_size=int(params['batchsize']), shuffle="batch", epochs=100,
                            validation_data=[[X_validpromoterSubseq, X_validhalflife], y_valid], callbacks=[earlystop_cb, check_cb])
        mse_history = result.history['val_mean_squared_error']
        mse = min(mse_history)
        best_file = os.path.join(params['datadir'], 'trainepoch.'+("%02d" % mse_history.index(mse))+'.h5')
        model = load_model(best_file)
        print 'Loaded results from:', best_file
        predictions_test = model.predict([X_testpromoterSubseq, X_testhalflife], batch_size=20).flatten()
        slope, intercept, r_value, p_value, std_err = stats.linregress(predictions_test, y_test)
        print 'Test R^2 = %.3f' % r_value**2

    return {'loss': mse, 'status': STATUS_OK }
    # except:
        # return {'loss': 1, 'status': STATUS_FAIL } # loss = 1 indicates a poor-performing model; reason model might fail include: incompatible parameters or insufficient memory resources available

if __name__ == '__main__':
    main()
