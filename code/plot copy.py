#%%
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import confusion_matrix
from helper_functions import train_test_split
from random import seed
import json
from sklearn.metrics import precision_recall_fscore_support
import matplotlib_venn as venn
from help_function import *
# %%
df_FY = pd.read_csv("/home/mddo/stage/M2S4/output/FY/diploid_/diploid_0/df/diploid_normal_KNN.csv")
df_CCD = pd.read_csv("/home/mddo/stage/M2S4/output/CCD/diploid_/diploid_0/df/normal_raw.csv")

FY_data = df_FY.loc[(df_FY["label"] == "non_ess")]
CCD_data = df_CCD.loc[(df_CCD["label"] == "non_ess")]

df_compare = pd.DataFrame()
df_compare["insertion_index_FY"] = FY_data["NI"]
df_compare["insertion_index_CCD"] = CCD_data["NI"]
sns.set(style="whitegrid")
# sns.boxplot(data=df_compare, showfliers = False)
sns.boxplot(x = "label", y = "NI", data = df_FY.loc[(df_FY["label"] == "ess")],showmeans=True, showfliers = False)

#%%
sns.boxplot(data=df_compare)
# %%
sns.boxplot(y = "label", x="hits_count", data=df_FY)
# %%
sns.boxplot(y = "label", x="hits_count", data=df_CCD)

# %%
sns.boxplot(x = "label", y="NI", data=df)

# %%
sns.boxplot(x = "predictions", y="NI", data=df)

# %%
sns.boxplot(x = "label", y="insertion_index", data=df)

# %%
sns.boxplot(x = "predictions", y="insertion_index", data=df)

# %%
sns.boxplot(x = "label", y="orf_len", data=df)

#%%

sns.boxplot(x = "predictions", y="orf_len", data=df)
# %%
df = pd.read_csv("/home/mddo/stage/M2S4/output/FY/diploid_/diploid_0/df/HFI_NI_PROM_new.csv")
sns.boxplot(x = "label", y="ratio_hits_prom", data=df)

# %%
df = pd.read_csv("/home/mddo/stage/M2S4/output/FY/diploid_/diploid_0/df/HFI_NI_PROM_KNN.csv")
sns.boxplot(x = "label", y="ratio_hits_prom", data=df)

# %%
df = pd.read_csv("/home/mddo/stage/M2S4/output/FY/diploid_/diploid_0/df/normal_KNN.csv")

# %%
data = df["label"].value_counts()

# %%
data.plot(kind="bar")

# %%
data_gene_count = df['label'].value_counts()
data_gene_count = pd.DataFrame(data_gene_count)
data_gene_count = data_gene_count.reset_index()
data_gene_count.columns = ["label","number of genes"]
print(data_gene_count)

sns.set(style="whitegrid")
#tips = sns.load_dataset("tips")
ax = sns.barplot(y = "number of genes", x = "label", data = data_gene_count, ci="sd")
for p in ax.patches:
    ax.annotate(str(p.get_height()), (p.get_x() * 1.005, p.get_height() * 1.005),size=20, 
                xytext=(180, 2),  # 3 points vertical offset
                textcoords="offset points", va='bottom')
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(18.5, 10.5)
plt.xlabel("Label",size = 20)
plt.ylabel("Number of gene",size = 20)
plt.xticks(size = 20)
plt.yticks(size = 20)
fig.savefig('test2png.png', dpi=100)

# %%
strain_names = ["FY","ABP","ACF","ACN","ACP","ADD","AND","APH","AVI","BBQ","BHH","BMK","CCD","CGQ","CHM","CIB","CLG","CNM","CNT","CPG","Sigma"]
acc_df_total = pd.DataFrame()
pre_df_total = pd.DataFrame()
for strain_name in strain_names:
    session_name = "test"
    type_data = "HFI_NI_PROM_NEW_removed"
    folder_number = 0

    accuracy_file = "/home/mddo/stage/M2S4/output/{}/accuracy/{}/accuracy_{}_{}.csv".format(strain_name, session_name, type_data,folder_number)
    # create_folder("/home/mddo/stage/M2S4/images/{}".format(strain_name))
    # create_folder("/home/mddo/stage/M2S4/images/{}/acc_and_precision/".format(strain_name))
    save_figure = "/home/mddo/stage/M2S4/images/{}/acc_and_precision/{}_{}.png".format(strain_name,type_data, folder_number)
    accuracy_df = pd.read_csv(accuracy_file)
    accuracy_df.columns = ["forest","accuracy","precision","recall","fscrore","total_tree"]

    accuracy = accuracy_df["accuracy"]
    precision = accuracy_df["precision"]
    total_tree = accuracy_df["total_tree"]
    acc_df_total["acc_{}".format(strain_name)] = accuracy_df["accuracy"]
    pre_df_total["prec_{}".format(strain_name)] = accuracy_df["precision"]


ax = plt.gca()
sns.set_style("whitegrid")

acc_df_total.plot(kind='line',ax=ax)
# accuracy_df.plot(kind='line',y='precision', color='red', ax=ax)
# accuracy_df = accuracy_df.drop(columns = ["total_tree"])
# accuracy_df.plot(kind="line")

fig = matplotlib.pyplot.gcf()
fig.set_size_inches(18.5, 10.5)
plt.rcParams["figure.figsize"] = (20,2)
plt.xticks(size = 14)
plt.yticks(size = 14)
plt.xlabel("Iteration", size = 18)
plt.ylabel("Accuracy", size = 18)
plt.ylim(0.85, 1)


# %%
strain_names = ["FY","ABP","ACF","ACN","ACP","ADD","AND","APH","AVI","BBQ","BHH","BMK","CCD","CGQ","CHM","CIB","CLG","CNM","CNT","CPG","Sigma"]
    # strain_names = ["FY"]
type_df = "HFI_NI_PROM_NEW_removed"
folder_number = 0
session_name = "test"
params = [strain_names, folder_number, session_name]
mean_score(type_df, params)

# %%
df = pd.read_csv("/home/mddo/stage/M2S4/output/FY/diploid_/diploid_0/df/HFI_NI_PROM_NEW_removed.csv")


# %%
sns.boxplot(x = df["NI"], y = "label", data = df)

# %%
df.loc[df["NI"] > 1.1, "NI"] = 1.1

# %%
sns.boxplot(x = df["NI"], y = "label", data = df)

# %%
# strain_names = ["FY","ABP","ACF","ACN","ACP","ADD","AND","APH","AVI","BBQ","BHH","BMK","CCD","CGQ","CHM","CIB","CLG","CNM","CNT","CPG","Sigma"]
# type_df = "HFI_NI_PROM_NEW_removed"
# folder_number = 0
# for strain_name in strain_names:
#     df_file = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/df/{}.csv".format(strain_name,folder_number,type_df)
#     df = pd.read_csv(df_file)
#     df.loc[df["NI"] > 1.1, "NI"] = 1.1
#     df.loc[df["hits_count"] >= 0.12, "hits_count"] = 0.12
#     df.loc[df["reads_count"] >= 0.04, "reads_count"] = 0.04
#     df.loc[df["insertion_index"] >= 0.025, "insertion_index"] = 0.025
#     df.loc[df["HFI_ratio"] >= 5, "HFI_ratio"] = 5
#     df.loc[df["NI_ratio"] >= 2, "NI_ratio"] = 2
#     df.loc[df["ratio_hits_prom"] >= 4, "ratio_hits_prom"] = 4
#     df.loc[df["reads_by_len"] >= 0.4, "reads_by_len"] = 0.4
#     df = df.drop(columns = ["orf_len"])
#     df.to_csv("/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/df/{}_1.csv".format(strain_name,folder_number,type_df), index=False)

# %%
df = pd.read_csv("/home/mddo/stage/M2S4/output/FY/accuracy/train/accuracy_HFI_NI_PROM_NEW_removed_1_0.csv")

# %%
df.mean()

# %%
