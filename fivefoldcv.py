import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import argparse

import xlsxwriter

from models import GraphConv, AE, LP
from utils import *
from sklearn.metrics import roc_curve, auc, roc_auc_score, precision_score, recall_score, f1_score, accuracy_score
parser = argparse.ArgumentParser()
parser.add_argument('--no-cuda', action='store_true', default=False,
                    help='Disables CUDA training.')
parser.add_argument('--seed', type=int, default=1, help='Random seed.')
parser.add_argument('--epochs', type=int, default=500,
                    help='Number of epochs to train.')
parser.add_argument('--lr', type=float, default=0.01,
                    help='Learning rate.')
parser.add_argument('--weight_decay', type=float, default=1e-5,
                    help='Weight decay (L2 loss on parameters).')
parser.add_argument('--hidden', type=int, default=256,
                    help='Dimension of representations')
parser.add_argument('--alpha', type=float, default=0.8,
                    help='Weight between lncRNA space and disease space')
parser.add_argument('--data', type=int, default=1, choices=[1,2],
                    help='Dataset')

args = parser.parse_args()
args.cuda = not args.no_cuda and torch.cuda.is_available()

set_seed(args.seed,args.cuda)
gdi, ldi, rnafeat, gl, gd = load_data(args.data,args.cuda)

class GNNq(nn.Module):
    def __init__(self):
        super(GNNq,self).__init__()
        self.gnnql = AE(rnafeat.shape[1],256,args.hidden)
        self.gnnqd = AE(gdi.shape[0],256,args.hidden)
    
    def forward(self,xl0,xd0):
        hl,stdl,xl = self.gnnql(gl,xl0)
        hd,stdd,xd = self.gnnqd(gd,xd0)
        return hl,stdl,xl,hd,stdd,xd

class GNNp(nn.Module):
    def __init__(self):
        super(GNNp,self).__init__()
        self.gnnpl = LP(args.hidden,ldi.shape[1])
        self.gnnpd = LP(args.hidden,ldi.shape[0])

    def forward(self,y0):
        yl,zl = self.gnnpl(gl,y0)
        yd,zd = self.gnnpd(gd,y0.t())
        return yl,zl,yd,zd

print("Dataset {}, 5-fold CV".format(args.data))

def criterion(output,target,msg,n_nodes,mu,logvar):
    if msg == 'disease':
        cost = F.binary_cross_entropy(output,target)
    else:
        cost = F.mse_loss(output,target)
    
    KL = -0.5 / n_nodes * torch.mean(torch.sum(
        1 + 2 * logvar - mu.pow(2) - logvar.exp().pow(2), 1))
    return cost + KL

def train(gnnq,gnnp,xl0,xd0,y0,epoch,alpha,i):
    losspl1 = []
    losspd1 = []
    lossp1 = []
    lossq1 = []
    beta0 = 1.0
    gamma0 = 1.0
    optp = torch.optim.Adam(gnnp.parameters(),lr=args.lr,weight_decay=args.weight_decay)
    optq = torch.optim.Adam(gnnq.parameters(),lr=args.lr,weight_decay=args.weight_decay)

    for e in range(epoch):

        gnnq.train()
        hl1, stdl1, xl1, hd1, stdd1, xd1 = gnnq(xl0, xd0)
        gnnq.train()
        hl2, stdl2, xl2, hd2, stdd2, xd2 = gnnq(xl0, xd0)
        gnnq.train()
        hl, stdl, xl, hd, stdd, xd = gnnq(xl2, xd2)
        lossql = criterion(xl, xl2,
                           "lncrna", gl.shape[0], hl, stdl)
        lossqd = criterion(xd, xd0,
                           "disease", gd.shape[0], hd, stdd)
        lossq = alpha * lossql + (1 - alpha) * lossqd + beta0 * e * F.mse_loss(
            torch.mm(hl, hd.t()), y0) / epoch
        optq.zero_grad()
        lossq1.append(lossq.item())
        lossq.backward()
        optq.step()
        gnnq.eval()
        with torch.no_grad():
            hl,_,_,hd,_,_ = gnnq(xl0,xd0)
        
        gnnp.train()
        yl,zl,yd,zd = gnnp(y0)
        losspl = F.binary_cross_entropy(yl,y0) + gamma0*e*F.mse_loss(zl,hl)/epoch
        losspd = F.binary_cross_entropy(yd,y0.t()) + gamma0*e*F.mse_loss(zd,hd)/epoch
        lossp = alpha*losspl + (1-alpha)*losspd
        losspl1.append(losspl.item())
        losspd1.append(losspd.item())
        lossp1.append(lossp.item())
        optp.zero_grad()
        lossp.backward()
        optp.step()


        with torch.no_grad():
            yl, _, yd, _ = gnnp(y0)
        if e%20 == 0:
            print('Epoch %d | Lossp: %.4f | Lossq: %.4f' % (e, lossp.item(),lossq.item()))
        r = pd.DataFrame(lossp1)
        r.to_csv('./output/lossp1{}.csv'.format(i))
        r1 = pd.DataFrame(lossq1)
        r.to_csv('./output/lossq1{}.csv'.format(i))
        r = pd.DataFrame(losspd1)
        r.to_csv('./output/losspd1{}.csv'.format(i))
        r = pd.DataFrame(losspl1)
        r.to_csv('./output/losspl1{}.csv'.format(i))
    return alpha*yl+(1-alpha)*yd.t()

def fivefoldcv(A,alpha=0.2):

    N = A.shape[0]
    idx = np.arange(N)
    np.random.shuffle(idx)
    res = torch.zeros(5,A.shape[0],A.shape[1])
    aurocl = np.zeros(5)
    auprl = np.zeros(5)
    for i in range(5):
        print("Fold {}".format(i+1))
        A0 = A.clone()
        k=i * N // 5
        k1=(i + 1) * N // 5
        for j in range(i*N//5,(i+1)*N//5):
            A1=A.shape[1]
            A0[idx[j],:] = torch.zeros(A.shape[1])
        
        gnnq = GNNq()
        gnnp = GNNp()
        if args.cuda:
            gnnq = gnnq.cuda()
            gnnp = gnnp.cuda()

        train(gnnq,gnnp,rnafeat,gdi.t(),A0,args.epochs,args.alpha,i)
        gnnq.eval()
        gnnp.eval()
        yli,_,ydi,_ = gnnp(A0)
        resi = alpha*yli + (1-alpha)*ydi.t()
        res[i] = resi
        
        if args.cuda:
            ymat = resi.cpu().detach().numpy()
        else:
            ymat = resi.detach().numpy()

        mean_tpr = 0.0
        mean_fpr = np.linspace(0, 1, 100)
        all_tpr = []
        all_fpr = []
        all_roc_auc = []
        tprs = []
        aucs = []

        global roc_auc_score
        wb = xlwt.Workbook(encoding='utf-8')
        results = []
        path = 'Dataset' + str(args.data)
        ldi = pd.read_csv(path + '/association_matrix.csv', header=None, encoding='gb18030').values
        y_true = ldi.flatten()
        ymat = ymat.flatten()
        fpr, tpr, rocth = roc_curve(y_true, ymat)

        auroc = auc(fpr, tpr)
        np.savetxt('./output/roc{}.txt'.format(i), np.vstack((fpr, tpr)), fmt='%10.5f', delimiter=',')
        precision, recall, prth = precision_recall_curve(y_true, ymat)
        aupr = auc(recall, precision)
        roc_sc = metrics.roc_auc_score(y_true, ymat)
        aupr_sc = metrics.average_precision_score(y_true, ymat)
        # apk_sc = rank_metrics.apk(actual, predicted, k=200)
        results.append([i, '%0.4f'%auroc, '%0.4f'%aupr, '%0.4f'%roc_sc, '%0.4f'%aupr_sc])
        r = pd.DataFrame(results)
        r.to_csv('./output/result{}.csv'.format(i))
        np.savetxt('./output/pr{}.txt'.format(i), np.vstack((recall, precision)), fmt='%10.5f', delimiter=',')
        print('AUROC= %.4f | AUPR= %.4f| roc_sc= %.4f| aupr_sc= %.4f' % (auroc, aupr, roc_sc, aupr_sc))
        rocdata = np.loadtxt('roc.txt', delimiter=',')
        prdata = np.loadtxt('pr.txt', delimiter=',')
        plt.figure()
        plt.plot(rocdata[0], rocdata[1])
        plt.plot(prdata[0], prdata[1])
        plt.show()

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
        aurocl[i] = auroc
        auprl[i] = aupr
        i=i+1

    mean_tpr /= 5
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(tprs, axis=0)
    plt.plot(mean_fpr, mean_tpr, color='b', label=r'Mean ROC (area=%0.2f)' % mean_auc, lw=2, alpha=.8)
    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    # plt.fill_between(mean_tpr, tprs_lower, tprs_upper, color='gray', alpha=.2)
    # plt.xlim([-0.05, 1.05])
    # plt.ylim([-0.05, 1.05])
    # plt.xlabel('False Positive Rate')
    # plt.ylabel('True Positive Rate')
    # plt.title('ROC')
    # plt.legend(loc='lower right')
    # plt.show()
    mean_tpr = mean_tpr.tolist()
    mean_fpr = mean_fpr.tolist()
    ws1 = xlsxwriter.Workbook('./output/ROC+' + '.xls')
    ws = ws1.add_worksheet('r1')
    ws.write(0, 0, 'Mean FPR')
    ws.write(1, 0, 'Mean TPR')
    for i in range(0, len(mean_fpr)):
        ws.write(0, i + 1, mean_fpr[i])
        ws.write(1, i + 1, mean_tpr[i])
    ws.write(2, 0, 'Mean ROC Area: %0.4f' % mean_auc)
    count = 3
    for num in range(0, len(all_tpr)):
        fold = num + 1
        ws.write(count, 0, 'FPR fold ' + str(fold))
        ws.write(count + 1, 0, 'TPR fold ' + str(fold))
        # break
        for i in range(0, len(all_tpr[num])):
            ws.write(count, i + 1, all_fpr[num][i])
            ws.write(count + 1, i + 1, all_tpr[num][i])
        ws.write(count + 2, 0, 'ROC Area: %0.4f' % all_roc_auc[num])
        count += 3
    print('OK!')
    ws1.close()
    # wb.save('./output/ROC+' + '.xls')
    ymat = res[auprl.argmax()]
    if args.cuda:
        return ymat.cpu().detach().numpy()
    else:
        return ymat.detach().numpy()

title = 'result--dataset'+str(args.data)
ymat = fivefoldcv(ldi,alpha=args.alpha)
title += '--fivefoldcv'
ymat = scaley(ymat)
np.savetxt('./output/'+title+'.txt',ymat,fmt='%10.5f',delimiter=',')
print("===Final result===")
show_auc(ymat,args.data)

