import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



import seaborn as sns
r =[0.896,0.9456,0.9455,0.9634]
#r = 2 * np.random.rand(100) #生成100个服从“0~1”均匀分布的随机样本值
theta = np.linspace( 2/5 * np.pi,4/5 * np.pi, 5)[:-1] #生成角度
area = 100 * np.array(r)**3 #面积
colors = ['darkorange','cornflowerblue','mediumslateblue', 'pink'] #颜色
ax = plt.subplot(111, projection='polar')
#projection为画图样式，除'polar'外还有'aitoff', 'hammer', 'lambert'等
c = ax.scatter(theta, r, c=colors, s=area,cmap='cool', alpha=0.75)
#ax.scatter为绘制散点图函数
plt.show()

import mpl_toolkits.mplot3d
# valueb=[0.9634,0.946,0.8533,0.7936,0.729,0.666,0.9007]
valueb=[0.9616,0.903,0.839,0.848,0.729,0.666,0.7996,0.926,0.882,0.819]
# valueb=[0.9173,0.8131,0.9353,0.9952,0.9784,0.9831,0.9088,0.8707,0.8857,0.8049,0.9443,0.8971,0.9568,0.8603,0.9935,0.9482,0.9769,0.8758,0.8548,0.8396,0.8178,0.8691,0.893,0.855,0.8599,0.828,0.9716,0.983,0.9684,0.9945,0.9263,0.9637,0.9217,0.9789,0.9306,0.9622,0.8105]
# valuesc=[0.87325,0.8025,0.8777,0.8675,0.8635,0.9516,0.9092,0.85706,0.8972,0.935,0.8988,0.9443,0.9527,0.7138,0.9475,0.973,0.6863,0.7029,0.6291,0.8106,0.7578,0.9544,0.9452,0.871,0.9581,0.8973,0.9245,0.8882,0.9465,0.9209,0.8282,0.9089,0.9482,0.9504,0.9524,0.911,0.9006]
# valuea = [0.9232,0.8233,0.9472,0.9952,0.981,0.9841,0.9088,0.8765,0.9236,0.8053,0.9443,0.8966,0.9618,
#        0.8618,0.9948,0.9518,0.9769,0.8758,0.8554,0.8426,0.8229,0.8751,0.9014,0.8649,0.8674,0.8347,0.9758,0.9879,0.9684,0.9945,0.9336,0.9666,0.9249,0.9797,0.9306,0.9713,0.8151]
#values = [0.09,-0.05,0.20,-0.02,0.08,0.09,0.03,0.027]
x = np.linspace(0, 1 * np.pi, 11)[:-1]
c = np.random.random(size=(10, 3))
# labels = ["AGO1","AGO2","AGO3","ALKBH5","AUF1","C17ORF85","C22ORF28","CAPRIN1","DGCR8","EIF4A3","EWSR1","FMRP"
#    ,"FOX2","FUS","FXR1","FXR2","HNRNPC","HUR","IGF2BP1","IGF2BP2","IGF2BP3","LIN28A","LIN28B"
#    ,"METTL3","MOV10","PTB","PUM2","QKI","SFRS1","TAF15","TDP43","TIA1","TIAL1","TNRC6","U2AF65","WTAP","ZC3HB"]
# labels=["GMNN2CD","ICFCDA","iCDA-CGR","KATZHCDA","GHICD","RWRHCD","CD-LNLP"]
labels=["GMNN2CD","GCNCDA","DWNN-RLS","KATZHCDA","GHICD","RWRHCD","CD-LNLP","NSL2CD","MRLDC","DeepDCR"]
fig = plt.figure(dpi=200)
plt.axes(polar=True)
# 获取当前的axes
print(plt.gca())
# 绘图
plt.bar(x, valueb, width=0.1, color=c, align='center')
plt.scatter(x, valueb, marker='o', c='black')
# 添加文本
plt.figtext(0.03, 0.7, s='', fontproperties='KaiTi', fontsize=22, rotation='vertical',
            verticalalignment='center', horizontalalignment='center')

plt.ylim(0.6, 0.95)


dataLength = 10
angles = np.linspace(0, 1 * np.pi, dataLength)
plt.thetagrids(angles * 180 / np.pi, labels, fontproperties='KaiTi', fontsize=10)
plt.grid(c='gray')
plt.show()
# tips = pd.read_csv(r'E:\博士论文相关\circRNA-disease\dp.csv')
# data=tips.values
# i=["GMNN2CD","ELM","SVM","RF","Recomm"]
# index=[chr(i) for i in range(52, 90)]
# df = pd.DataFrame(data,
#                   index=i,
#                   columns=["Dataset1","Dataset2","Dataset3","Dataset4","Dataset5"])
# # plt.clf()
# plt.figure(dpi=200)
#
# sns.heatmap(data=df,  cmap=plt.get_cmap('Greens'),annot=True,fmt=".4f",
#             annot_kws={'size':7,'weight':'normal', 'color':'blue'})
#
# #plt.title("使用matplotlib中的颜色盘：cmap=plt.get_cmap('Greens')")
# plt.show()
# 数据格式如下
# disease = ["Meningitis, Pneumococcal", "Meningitis, Pneumococcal", "Mycetoma", "Botulism", "Botulism2", "Botulism3"]
# id = ["C01.252.200.500.600", "C08.345.654.570", "C01.252.410.040.692.606", "C01.252.410.222.151", "C01.252.410.222.151", "C03.252.410.222.151"]
# losp1 = pd.read_csv('./epoch/lossp1.csv', header=0, encoding='gb18030').values
# losq1=pd.read_csv('./epoch/lossq1.csv', header=0, encoding='gb18030').values
# losp2 = pd.read_csv('./epoch/lossp2.csv', header=0, encoding='gb18030').values
# losq2=pd.read_csv('./epoch/lossq2.csv', header=0, encoding='gb18030').values
# losp3 = pd.read_csv('./epoch/lossp3.csv', header=0, encoding='gb18030').values
# losq3=pd.read_csv('./epoch/lossq3.csv', header=0, encoding='gb18030').values
# losp4 = pd.read_csv('./epoch/lossp4.csv', header=0, encoding='gb18030').values
# losq4=pd.read_csv('./epoch/lossq4.csv', header=0, encoding='gb18030').values
# losp5= pd.read_csv('./epoch/lossp5.csv', header=0, encoding='gb18030').values
# losq5=pd.read_csv('./epoch/lossq5.csv', header=0, encoding='gb18030').values
#
# # rocdata1 = np.loadtxt('./epoch/roc1.txt', delimiter=',')
# # prdata1 = np.loadtxt('./epoch/pr1.txt', delimiter=',')
# # rocdata2 = np.loadtxt('./epoch/roc2.txt', delimiter=',')
# # prdata2 = np.loadtxt('./epoch/pr2.txt', delimiter=',')
# # rocdata3 = np.loadtxt('./epoch/roc3.txt', delimiter=',')
# # prdata3 = np.loadtxt('./epoch/pr3.txt', delimiter=',')
# # rocdata4 = np.loadtxt('./epoch/roc4.txt', delimiter=',')
# # prdata4= np.loadtxt('./epoch/pr4.txt', delimiter=',')
# # rocdata = np.loadtxt('./epoch/roc5.txt', delimiter=',')
# # prdata = np.loadtxt('./epoch/pr5.txt', delimiter=',')
# plt.figure(dpi=200)
# plt.plot(losq1[:,0], losq1[:,1],label='Dataset1-q(Inference)')
# plt.plot(losq1[:,0], losp1[:,1],label='Dataset2-p(Learning)')
# plt.plot(losq2[:,0], losq2[:,1],label='Dataset2-q(Inference)')
# plt.plot(losq2[:,0], losp2[:,1],label='Dataset2-p(Learning)')
# plt.plot(losq1[:,0], losq3[:,1],label='Dataset3-q(Inference)')
# plt.plot(losq1[:,0], losp3[:,1],label='Dataset3-p(Learning)')
# plt.plot(losq2[:,0], losq4[:,1],label='Dataset4-q(Inference)')
# plt.plot(losq2[:,0], losp4[:,1],label='Dataset4-p(Learning)')
# plt.plot(losq2[:,0], losq5[:,1],label='Dataset5-q(Inference)')
# plt.plot(losq2[:,0], losp5[:,1],label='Dataset5-p(Learning)')
# # plt.plot(rocdata1[0], rocdata1[1],linewidth=2,label='Dataset1')
# # plt.plot(rocdata2[0], rocdata2[1],linewidth=2,label='Dataset2')
# # plt.plot(rocdata3[0], rocdata3[1],linewidth=2,label='Dataset3')
# # plt.plot(rocdata4[0], rocdata4[1],linewidth=2,label='Dataset4')
# # plt.plot(rocdata[0], rocdata[1],linewidth=2,label='Dataset')
#
# # plt.plot(prdata1[0], prdata1[1],linewidth=2,label='Dataset1')
# # plt.plot(prdata2[0], prdata2[1],linewidth=2,label='Dataset2')
# # plt.plot(prdata3[0], prdata3[1],linewidth=2,label='Dataset3')
# # plt.plot(prdata4[0], prdata4[1],linewidth=2,label='Dataset4')
# # plt.plot(prdata[0], prdata[1],linewidth=2,label='Dataset5')
# plt.legend()  # 让图例生效
# # plt.xlabel("Recall") #X轴标签
# # plt.ylabel("Precission") #Y轴标签
# plt.xlabel("epoch") #X轴标签
# plt.ylabel("loss") #Y轴标签
# plt.show()

# print("开始读取数据")
# # 读取数据
# meshid = pd.read_csv('data/MeSHID.csv', header=0)
# disease = meshid['disease'].tolist()
# id = meshid['ID'].tolist()
#
# meshdis = pd.read_csv('data/Mesh_disease.csv', header=0)
# unique_disease = meshdis['C1'].tolist()
#
# # 初始化字典，有重复也没关系
# for i in range(len(disease)):
#     disease[i] = {}
#
# print("开始计算每个病的DV")
# # 计算每个病的DV，又重复也没关系，之后再合并
#
# for i in range(len(disease)):
#     l=len(id[i])
#     if len(id[i]) > 3:
#         disease[i][id[i]] = 1
#         id[i] = id[i][:-4]
#         # print(disease[i])
#         if len(id[i]) > 3:
#             disease[i][id[i]] = round(1 * 0.8, 5)
#             id[i] = id[i][:-4]
#             # print(disease[i])
#             if len(id[i]) > 3:
#                 disease[i][id[i]] = round(1 * 0.8 * 0.8, 5)
#                 id[i] = id[i][:-4]
#                 # print(disease[i])
#                 if len(id[i]) > 3:
#                     disease[i][id[i]] = round(1 * 0.8 * 0.8 * 0.8, 5)
#                     id[i] = id[i][:-4]
#                     # print(disease[i])
#                     if len(id[i]) > 3:
#                         disease[i][id[i]] = round(1 * 0.8 * 0.8 * 0.8 * 0.8, 5)
#                         id[i] = id[i][:-4]
#                         # print(disease[i])
#                         if len(id[i]) > 3:
#                             disease[i][id[i]] = round(1 * 0.8 * 0.8 * 0.8 * 0.8 * 0.8, 5)
#                             id[i] = id[i][:-4]
#                             # print(disease[i])
#                             if len(id[i]) > 3:
#                                 disease[i][id[i]] = round(1 * 0.8 * 0.8 * 0.8 * 0.8 * 0.8 * 0.8, 5)
#                                 id[i] = id[i][:-4]
#                                 # print(disease[i])
#                                 if len(id[i]) > 3:
#                                     disease[i][id[i]] = round(1 * 0.8 * 0.8 * 0.8 * 0.8 * 0.8 * 0.8 * 0.8, 5)
#                                     id[i] = id[i][:-4]
#                                     # print(disease[i])
#                                 else:
#                                     disease[i][id[i][:3]] = round(1 * 0.8 * 0.8 * 0.8 * 0.8 * 0.8 * 0.8 * 0.8, 5)
#                                     # print(disease[i])
#                             else:
#                                 disease[i][id[i][:3]] = round(1 * 0.8 * 0.8 * 0.8 * 0.8 * 0.8 * 0.8, 5)
#                                 # print(disease[i])
#                         else:
#                             disease[i][id[i][:3]] = round(1 * 0.8 * 0.8 * 0.8 * 0.8 * 0.8, 5)
#                             # print(disease[i])
#                     else:
#                         disease[i][id[i][:3]] = round(1 * 0.8 * 0.8 * 0.8 * 0.8, 5)
#                         # print(disease[i])
#                 else:
#                     disease[i][id[i][:3]] = round(1 * 0.8 * 0.8 * 0.8, 5)
#                     # print(disease[i])
#             else:
#                 disease[i][id[i][:3]] = round(1 * 0.8 * 0.8, 5)
#                 # print(disease[i])
#         else:
#             disease[i][id[i][:3]] = round(1 * 0.8, 5)
#             # print(disease[i])
#     else:
#         disease[i][id[i][:3]] = 1
#         # print(disease[i])
#
# print("合并相同的病不同ID的DV")
#
# # 合并相同的病不同ID的DV
#
# unique_disease = meshdis['C1'].tolist()
#
# # 这个name用来判断
# disease_name = meshid['disease'].tolist()
# unique_disease_name = meshdis['C1'].tolist()
#
# for i in range(len(unique_disease)):
#     unique_disease[i] = {}
#     for j in range(len(disease_name)):
#         if unique_disease_name[i] == disease_name[j]:
#             unique_disease[i].update(disease[j])
#
#
# similarity = np.zeros([len(unique_disease_name), len(unique_disease_name)])
#
# # print(similarity)
#
# print("计算相似度")
#
# for m in range(len(unique_disease_name)):
#     for n in range(len(unique_disease_name)):
#         u=unique_disease[m].values()
#         u1=unique_disease[n].values()
#         denominator = sum(unique_disease[m].values()) + sum(unique_disease[n].values())
#         numerator = 0
#         u2=unique_disease[m].items()
#         u3=unique_disease[n].keys()
#         for k, v in unique_disease[m].items():
#             if k in unique_disease[n].keys():
#                 numerator += v + unique_disease[n].get(k)
#         similarity[m, n] = round(numerator/denominator, 5)
#
# # print(similarity)
# print("保存结果")
#
# # 保存结果
#
# result = pd.DataFrame(similarity)
# result.to_csv('output/similarity1.csv')
