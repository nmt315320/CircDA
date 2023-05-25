import numpy as np
import torch
import matplotlib.pyplot as plt
import argparse
from sklearn import metrics
#from skchem.metrics import bedroc_score
from sklearn.preprocessing import minmax_scale,scale
import xlwt
import math
from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_curve,roc_auc_score,average_precision_score,precision_recall_curve,auc
import pandas as pd
from scipy import interp
from sklearn.metrics import roc_curve, auc, roc_auc_score, precision_score, recall_score, f1_score, accuracy_score
# from skchem.metrics import bedroc_score
def scaley(ymat):
    return (ymat-ymat.min())/ymat.max()

def set_seed(seed,cuda):
    np.random.seed(seed)
    torch.manual_seed(seed)
    if cuda:
        torch.cuda.manual_seed(seed)

def load_data(data,cuda):
    path = 'Dataset'+str(data)
    gdi = pd.read_csv(path + '/disease_GaussianSimilarity.csv',header=None, encoding='gb18030').values
    ldi = pd.read_csv(path + '/association_matrix.csv',header=None, encoding='gb18030').values
    rnafeat =pd.read_csv(path + '/rna_GaussianSimilarity.csv',header=None, encoding='gb18030').values
    rnafeat = minmax_scale(rnafeat,axis=0)
    gdit = torch.from_numpy(gdi).float()
    ldit = torch.from_numpy(ldi).float()
    rnafeatorch = torch.from_numpy(rnafeat).float()
    gl = norm_adj(rnafeat)
    gd = norm_adj(gdi.T)
    if cuda:
        gdit = gdit.cuda()
        ldit = ldit.cuda()
        rnafeatorch = rnafeatorch.cuda()
        gl = gl.cuda()
        gd = gd.cuda()
    
    return gdit, ldit, rnafeatorch, gl, gd

def neighborhood(feat,k):
    # compute C
    featprod = np.dot(feat.T,feat)
    smat = np.tile(np.diag(featprod),(feat.shape[1],1))
    dmat = smat + smat.T - 2*featprod
    dsort = np.argsort(dmat)[:,1:k+1]
    C = np.zeros((feat.shape[1],feat.shape[1]))
    for i in range(feat.shape[1]):
        for j in dsort[i]:
            C[i,j] = 1.0
    
    return C

def normalized(wmat):
    deg = np.diag(np.sum(wmat,axis=0))
    degpow = np.power(deg,-0.5)
    degpow[np.isinf(degpow)] = 0
    W = np.dot(np.dot(degpow,wmat),degpow)
    return W

def norm_adj(feat):
    C = neighborhood(feat.T,k=10)
    norm_adj = normalized(C.T*C+np.eye(C.shape[0]))
    g = torch.from_numpy(norm_adj).float()
    return g

def show_auc(ymat,data):
    mean_tpr = 0.0
    mean_fpr = np.linspace(0, 1, 100)
    all_tpr = []
    all_fpr = []
    all_roc_auc = []
    tprs=[]
    aucs = []

    global roc_auc_score
    wb = xlwt.Workbook(encoding='utf-8')
    results = []
    path = 'Dataset'+str(data)
    ldi = pd.read_csv(path + '/association_matrix.csv', header=None, encoding='gb18030').values
    y_true = ldi.flatten()
    ymat = ymat.flatten()
    fpr,tpr,rocth = roc_curve(y_true,ymat)

    auroc = auc(fpr,tpr)
    np.savetxt('roc.txt',np.vstack((fpr,tpr)),fmt='%10.5f',delimiter=',')
    precision,recall,prth = precision_recall_curve(y_true,ymat)
    aupr = auc(recall,precision)
    roc_sc = metrics.roc_auc_score(y_true,ymat)
    aupr_sc = metrics.average_precision_score(y_true,ymat)
    #apk_sc = rank_metrics.apk(actual, predicted, k=200)

    np.savetxt('pr.txt',np.vstack((recall,precision)),fmt='%10.5f',delimiter=',')
    print('AUROC= %.4f | AUPR= %.4f| roc_sc= %.4f| aupr_sc= %.4f'% (auroc,aupr,roc_sc,aupr_sc))
    # rocdata = np.loadtxt('roc.txt',delimiter=',')
    # prdata = np.loadtxt('pr.txt',delimiter=',')
    # plt.figure()
    # plt.plot(rocdata[0],rocdata[1])
    # plt.plot(prdata[0],prdata[1])
    # # plt.show()

    aucc = roc_auc_score(y_true, ymat)
    aucs.append(aucc)

    fpr, tpr, thresholds = roc_curve(y_true, ymat)
    tprs.append(interp(mean_fpr, fpr, tpr))
    mean_tpr += interp(mean_fpr, fpr, tpr)
    mean_tpr[0] = 0.0
    auc(fpr, tpr)
    roc_auc = auc(fpr, tpr)
    all_tpr.append(tpr.tolist())
    all_fpr.append(fpr.tolist())
    all_roc_auc.append(roc_auc)
    # print('AUC: ', aucc)
    # results.append([i,  '%0.4f' % recall,
    #                   '%0.4f' % aucc,'%0.4f' % precision])
    # plt.plot(fpr, tpr, lw=1, alpha=0.3, label='ROC fold %d(area=%0.2f)' % (i, roc_auc))
    # r = pd.DataFrame(results)
    # r.to_csv('result{}.csv'.format(i))

    return auroc,aupr