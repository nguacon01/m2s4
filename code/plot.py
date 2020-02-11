#%%
import numpy as numpy
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import confusion_matrix
from helper_functions import train_test_split

#%%
df = pd.read_csv("/home/mddo/stage/M2S4/output/FY/diploid/diploid_0/dataframe.csv")
df.shape
train_df, test_df = train_test_split(df, 0.2)
# data = test_df['label'].value_counts()
# print(data)
# data.plot(kind="bar",color = ["blue","orange"])
# plt.xticks(rotation=0,size=16)
# plt.yticks(size=14)
# plt.ylabel("numbers of gene",size=14)
# plt.show()
# ax = sns.barplot(df['label'].value_counts())



# ax = sns.scatterplot(x="HFI_ratio", y="NI_ratio", hue="label", data=df)

# sns.heatmap([train_df.hits_count,train_df.read_count],annot=True,fmt="d")

# %%
sns.scatterplot(y='HFI_ratio', x='NI_ratio', data = df, hue='label')

# # %%
# df = pd.read_csv("/home/mddo/stage/M2S4/output/output_predictions/output_predictions_epoche_0.csv")
# confusion_matrix = pd.crosstab(df['label'],df['predictions'], rownames = ['Actual'], colnames=['Predict'])
# sns.heatmap(confusion_matrix,annot=True,annot_kws={"size": 16,"color":"white"})
# plt.show()



# %%
# df['hits_count'].plot(kind='bar')

# # %%
# ver2_df = pd.read_csv("/home/mddo/stage/M2S4/output/FY/dataframe.csv")

# print(ver2_df.head(10))# %%


# # %%


# # %%
# sns.scatterplot(y='insertion_index', x='ratio_hits_prom', data = ver2_df, hue='label')
# plt.yticks(size = 14)
# plt.xticks(size = 14)
# plt.xlabel("ratio hits prom", size = 16)
# plt.ylabel("hits count in 100 bp prom", size = 16)
# plt.rcParams["figure.figsize"] = (30,4)

# %%
dataFrame = pd.read_csv("/home/mddo/stage/M2S4/logs/report.out")
# %%
# %%
sns.lineplot(data=pd.DataFrame(dataFrame["accuracy_HFI_NI"]), palette="tab10", linewidth=2.5, color = "orange",markers="o")

# plt.rcParams["figure.figsize"] = (30,4)
# sns.set(style="whitegrid")
# dataFrame.plot(kind="line")
plt.rcParams["figure.figsize"] = (30,4)
sns.set(style="whitegrid")

# %%

df_HFI_NI = pd.read_csv("/home/mddo/stage/M2S4/output/FY/accuracy/full_5k/50/training/HFI_NI/accuracy_HFI_NI.out")
df_HFI_NI.columns = ["n_tree", "n_feature", "n_max_depth", "n_bootstrap", "accuracy"]
df_HFI_NI_PROM = pd.read_csv("/home/mddo/stage/M2S4/output/FY/accuracy/full_5k/50/training/HFI_NI_PROM/accuracy_HFI_NI_PROM.out")
df_HFI_NI_PROM.columns = ["n_tree", "n_feature", "n_max_depth", "n_bootstrap", "accuracy"]
df_normal = pd.read_csv("/home/mddo/stage/M2S4/output/FY/accuracy/full_5k/50/training/normal/accuracy_normal.out")
df_normal.columns = ["n_tree", "n_feature", "n_max_depth", "n_bootstrap", "accuracy"]

# %%
df_HFI_NI_acc = pd.DataFrame(df_HFI_NI["accuracy"])

# %%
df_HFI_NI_PROM_acc = pd.DataFrame(df_HFI_NI_PROM["accuracy"])

# %%
df_normal_acc = pd.DataFrame(df_normal["accuracy"])

#%%
merge_df = pd.concat([df_HFI_NI_acc,df_HFI_NI_PROM_acc,df_normal_acc],axis=1, join='inner')
merge_df.columns = ["acc_HFI_NI","acc_HFI_NI_PROM","acc_normal"]
# %%
print(merge_df)

sns.set(style="whitegrid")
merge_df.plot(kind="line")
plt.rcParams["figure.figsize"] = (20,5)
plt.ylabel("accuracy")
plt.xlabel("epoches")
# %%


# %%
