# -*- coding: utf-8 -*-
"""
Created on Thu Jul  4 21:17:19 2019

@author: ZQQ
"""
import numpy as np
from sklearn.ensemble import AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import LeaveOneOut
from sklearn import datasets

iris = datasets.load_iris()
X = iris.data
y = iris.target

##变为2分类
X, y = X[y != 2], y[y != 2]
# 这个地方可以加上上一篇博客的随机打乱数据操作
loo = LeaveOneOut()
loo.get_n_splits(X)
print("交叉验证次数：", loo.get_n_splits(X))  # 输出为100，--->进行100折,也就是留一

y_pred = []
for train_index, test_index in loo.split(X):
    # print("train:", train_index, "TEST:", test_index) # 索引
    X_train, X_test = X[train_index], X[test_index]
    y_train, y_test = y[train_index], y[test_index]
    # print(X_train, X_test, y_train, y_test)

    # 调用、训练模型
    model_bdt = AdaBoostClassifier(DecisionTreeClassifier(max_depth=2), algorithm="SAMME", n_estimators=10)
    model_bdt.fit(X_train, y_train)

    # 预测
    x_test_pred = model_bdt.predict(X_test)

    # print(x_test_pred)
    y_pred.append(x_test_pred)  # 当前预测值添加到列表中

from sklearn.metrics import roc_curve, auc

y_pred = np.array(y_pred)  # list to array
fpr, tpr, threshold = roc_curve(y, y_pred)  # 计算真正率和假正率
roc_auc = auc(fpr, tpr)  # 计算auc的值

import matplotlib.pyplot as plt

lw = 2  # 定义线条宽度
plt.figure(figsize=(8, 5))
plt.plot(fpr, tpr, color='darkorange',
         lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)  ###假正率为横坐标，真正率为纵坐标做曲线
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic example')
plt.legend(loc="lower right")
plt.savefig('loocv.png', dpi=600)  # 以600大批保存图片
plt.show()
