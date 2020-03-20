import pickle, sys

file = sys.argv[1]
trials = pickle.load(open(file, "rb"))
alldicts = trials.trials
loss = trials.losses()

for i in range(len(loss)):
    print(str(int(alldicts[i]['misc']['vals']['leftpos'][0]))+'\t'+str(int(alldicts[i]['misc']['vals']['rightpos'][0]))+'\t'+str(loss[i])+'\t'+str(alldicts[i]['misc']['vals']))
