# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 08:43:50 2019

@author: caseyp
"""

# import libraries
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split, GridSearchCV, cross_val_score, ShuffleSplit
from sklearn.metrics import mean_squared_error, r2_score
from scipy.io import loadmat, savemat
from scipy.stats import zscore

# load data
homeDir='/home/cpaq3/Desktop/DMN/'
mat_file = loadmat(homeDir + '/output/bigbrain_features.mat')

# split data
X_train,X_test,y_train,y_test = train_test_split(zscore(mat_file['X']), mat_file['y'], test_size=0.3)
y_train = y_train.ravel()
y_test = y_test.ravel()

# Perform Grid-Search
gsc = GridSearchCV(
    estimator=RandomForestRegressor(),
    param_grid={
             'max_depth': range(3,10),
            'n_estimators': range(5, 40, 5),
            },
    cv=5, scoring='neg_mean_squared_error', verbose=0,n_jobs=-1)
grid_result = gsc.fit(X_train, y_train)
best_params = grid_result.best_params_
    
# build full model
rfr = RandomForestRegressor(max_depth=best_params["max_depth"], 
                            n_estimators=best_params["n_estimators"],
                            random_state=False, verbose=False)

# feature selection
scores = np.zeros([10,X.shape[1]])
acc_full = np.zeros(10)
rs = ShuffleSplit(n_splits=5, test_size=0.3, random_state=0)
x = 0
for train_idx, test_idx in rs.split(X_train):
    X_train2, X_test2 = X_train[train_idx], X_train[test_idx]
    y_train2, y_test2 = y_train[train_idx], y_train[test_idx]
    rfr.fit(X_train2, y_train2)
    acc_full[x] = mean_squared_error(y_test2, rfr.predict(X_test2))
    for i in range(0,X.shape[1]):
        X_t = X_test2.copy()
        np.random.shuffle(X_t[:, i])
        shuff_acc = mean_squared_error(y_test2, rfr.predict(X_t))
        tmp=(acc_full[x]-shuff_acc)/acc_full[x]
        scores[x,i]=tmp
    x = x + 1

# select features that when shuffled are associated with greater than average decrease in variance explained
sel =  np.mean(scores,0)<np.mean(scores) 

# retrain and test model with selected features
rfr_slim = RandomForestRegressor(max_depth=best_params["max_depth"], 
                            n_estimators=best_params["n_estimators"],
                            random_state=False, verbose=False)              
rfr_slim.fit(X_train[:,sel], y_train)
acc_sel = r2_score(y_test, rfr_slim.predict(X_test[:,sel]))
mse_sel = mean_squared_error(y_test, rfr_slim.predict(X_test[:,sel]))

np.savetxt(homeDir + '/output/random_forest.csv', scores.values(), delimiter=",")