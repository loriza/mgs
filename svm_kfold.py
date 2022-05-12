# -*- coding: utf-8 -*-
"""
Created on Wed Nov 24 18:05:18 2021

@author: lymon
"""

import numpy as np
import pandas as pd
import argparse
from scipy import interp
import matplotlib.pyplot as plt

from sklearn import svm
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold


parser = argparse.ArgumentParser()
parser.add_argument("-i", help="-i feature table, required", required=True)
parser.add_argument("-l", help="-l labels")
parser.add_argument("-o", help="-o prefix of outfile, required", required=True)
parser.add_argument(
    "-k",type = int, help="-k k-fold , default 10", default="10")

args = parser.parse_args()

# #############################################################################
# Data IO and generation

def load_data(InputPath):
    df = pd.read_table(InputPath, sep='\t',header = 0, index_col=0)
    #print(df)
    input = np.array(df)

    return input

 #############################################################################
# Classification and ROC analysis
C=1.0 #惩罚参数，pred val越高泛化能力越弱
Kernel='rbf'
degree=3
gamma='auto'
coef0=0.0
shrinking=True
tol=0.001
cache_size=200
class_weight=None
verbose=False
max_iter=-1
#decision_function_shape='ovr'
RandomState=123


         
def cv_ROC(featurePath,labelPath,output,cvfold):  
    
    X, y = load_data(featurePath), load_data(labelPath)

    cv = StratifiedKFold(n_splits=cvfold)
    sv = svm.SVC(kernel=Kernel, probability=True,random_state=RandomState)  
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)
    
    i = 0
    for train, test in cv.split(X, y):
        # returns the false positive rates for each threshold, true positive rates for each threshold and thresholds.
        probas = sv.fit(X[train], y[train].ravel()).predict_proba(X[test])     
        
        # Compute ROC curve and area the curve
        fpr, tpr, thresholds = roc_curve(y[test].ravel(), probas[:, 1],pos_label =1)
        print(fpr)
        tprs.append(np.interp(mean_fpr, fpr, tpr))
        tprs[-1][0] = 0.0
        roc_auc = auc(fpr, tpr)
        aucs.append(roc_auc)
        plt.plot(fpr, tpr, lw=1, alpha=0.3,
                 label='ROC fold %d (AUC = %0.2f)' % (i+1, roc_auc))
    
        i = i + 1

    plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
              alpha=.8) 
    #compute the average ROC curve and AUC
    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    plt.plot(mean_fpr, mean_tpr, color='b',
             label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
             lw=2, alpha=.8)
    
    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)

    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc="lower right")
   
    plt.savefig(output+"_svm_kfold_roc.png",bbox_inches ='tight')
  
    return plt




if __name__ == '__main__':
    # Run classifier with cross-validation 
    
    cv_ROC(args.i,args.l,args.o,args.k)
