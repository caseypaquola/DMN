# import libraries
import numpy as np
from sklearn.multiclass import OneVsOneClassifier
from sklearn.svm import SVC
from sklearn import svm, datasets
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
import scipy.io as io
from scipy import interp

# load data
projDir='/data_/mica1/03_projects/casey/sandbox1/9_DMN/'
mat = io.loadmat(projDir + "output/bigbrain_features_subdmn.mat")
X = mat["Xcv_eq"]
ycv = mat["ycv_eq"]
n_classes = 3
folds = X.shape[2]
obs = ycv.shape[0]

# predict in each fold
y_pred = np.zeros([round(obs/4),folds])
y_comp = np.zeros([round(obs/4), folds])
fpr = np.zeros([3,folds,3])
tpr = np.zeros([3,folds,3])
for i in range(0,folds):
    X_train, X_test, y_train, y_test = train_test_split(X[:,:,i],
        ycv[:,i], test_size=.25,
        random_state=i)
    classifier = OneVsOneClassifier(svm.SVC(kernel='linear', probability=True))
    classifier.fit(X_train, y_train)
    y_pred[:,i] = classifier.predict(X_test)
    y_comp[:,i] = y_test
    y_pred_bin = label_binarize(y_pred[:,i], classes=[1, 2, 3])
    y_test_bin = label_binarize(y_test, classes=[1, 2, 3])
    for j in range(n_classes):
        fpr[:,i,j], tpr[:,i,j], _ = roc_curve(y_test_bin[:, j], y_pred_bin[:, j])

io.savemat(projDir + "output/bigbrain_subdmn_pred.mat", 
    {"y_pred": y_pred, "y_comp": y_comp, "fpr": fpr, "tpr": tpr})

# predict in each fold, and do leave-one-out feature
y_pred = np.zeros([round(obs/4),folds, X.shape[1]]) # note: only using the gradients
y_comp = np.zeros([round(obs/4),folds])
for i in range(0,folds):
    X_train, X_test, y_train, y_test = train_test_split(X[:,:,i], 
        ycv[:,i], test_size=.25,
        random_state=0)
    y_comp[:,i] = y_test
    for j in range(0, X.shape[1]):
        classifier = OneVsOneClassifier(svm.SVC(kernel='linear', probability=True))
        X_train_loo = np.delete(X_train, j, 1)
        classifier.fit(X_train_loo, y_train)
        X_test_loo= np.delete(X_test, j, 1)
        y_pred[:,i,j] = classifier.predict(X_test_loo)

io.savemat(projDir + "output/bigbrain_subdmn_pred_loo.mat", {"y_pred": y_pred, "y_comp": y_comp})        