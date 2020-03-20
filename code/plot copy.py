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
# %%
df = pd.read_csv("/home/mddo/stage/M2S4/output/FY/diploid_/diploid_0/df/normal_KNN.csv")

# %%
sns.boxplot(x = "label", y="hits_count", data=df)
# %%
sns.boxplot(x = "label", y="reads_count", data=df)

# %%
sns.boxplot(x = "label", y="NI", data=df)

# %%
sns.boxplot(x = "label", y="HFI", data=df)

# %%
sns.boxplot(x = "label", y="insertion_index", data=df)

# %%
sns.boxplot(x = "label", y="hits_count_pro", data=df)

# %%
sns.boxplot(x = "label", y="orf_len", data=df)

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
strain_names = ["FY"]
acc_df_total = pd.DataFrame()
pre_df_total = pd.DataFrame()
for strain_name in strain_names:
    session_name = "test"
    type_data = "normal_KNN_removed"
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
    acc_df_total["prec_{}".format(strain_name)] = accuracy_df["precision"]
    print(acc_df_total)

    print(pre_df_total)

ax = plt.gca()
sns.set_style("whitegrid")

acc_df_total.plot(kind='line',ax=ax)
# accuracy_df.plot(kind='line',y='precision', color='red', ax=ax)
# accuracy_df = accuracy_df.drop(columns = ["total_tree"])
# accuracy_df.plot(kind="line")

fig = matplotlib.pyplot.gcf()
fig.set_size_inches(18.5, 5.5)
plt.rcParams["figure.figsize"] = (20,2)
# plt.title("{} prediction accuracy and precision".format(strain_name), size = 20)
plt.xticks(size = 14)
plt.yticks(size = 14)
plt.xlabel("Iteration", size = 18)
plt.ylabel("Accuracy", size = 18)


# %%
