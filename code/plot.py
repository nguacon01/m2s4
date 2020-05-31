#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import os
import seaborn as sns
#%%
strains_name = ["FY","ABP","ACF","ACN","ACP","ADD","AND","APH","AVI","BBQ","BHH","BMK","CCD","CGQ","CHM","CIB","CLG","CNM","CNT","CPG","Sigma"]

df = pd.read_csv("/home/mddo/stage/M2S4/output/FY/diploid_/diploid_0/df/HFI_NI_PROM_NEW.csv_removed.csv")


df_ess = df.loc[df["label"] == "ess"]
print(df_ess.shape)
df_noness = df.loc[df["label"] == "non_ess"]
print(df_noness.index)
# bootstrap_indices = np.random.choice(len(df_noness),1838, replace=True)
# print(len(bootstrap_indices))
# for strain_name in strains_name:
#     df_ = pd.read_csv("/home/mddo/stage/M2S4/output/{}/diploid_/diploid_0/df/HFI_NI_PROM_NEW.csv_removed.csv".format(strain_name))
#     print("/home/mddo/stage/M2S4/output/{}/diploid_/diploid_0/df/HFI_NI_PROM_NEW.csv_removed.csv".format(strain_name))
#     # bootstrap_indices = np.random.randint(low=0,high=len(train_df), size=n_bootstrap)
#     df_bootstrapped = df_.iloc[bootstrap_indices]
#     print(df_bootstrapped.head(10))

#     df_concat = pd.concat([df_ess,df_bootstrapped])
#     print(df_concat.shape)
#     df_concat.to_csv("/home/mddo/stage/M2S4/output/{}/diploid_/diploid_0/df/balance_HFI_NI_PROM_NEW.csv".format(strain_name), index=False)
#     label_df = pd.DataFrame()
#     label_df["orf"] = df_concat["orf"]
#     label_df["label"] = df_concat["label"]
#     label_df.to_csv("label_balance.csv", index=False)
#     print(label_df.shape)


# %%


# %%
df = pd.read_csv("/home/mddo/stage/M2S4/output/FY/diploid_/diploid_0/df/balance_HFI_NI_PROM_NEW.csv")

# %%
value_count = df["label"].value_counts()

# %%
train_data_gene_count = pd.DataFrame(value_count)
train_data_gene_count = train_data_gene_count.reset_index()
train_data_gene_count.columns = ["label","number of gene"]
print(train_data_gene_count)


sns.set(style="whitegrid")
#tips = sns.load_dataset("tips")
ax = sns.barplot(y = "number of gene", x = "label", data = train_data_gene_count, ci="sd")
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


# %%
strain_names = ["FY"]
# strain_names = ["FY"]
acc_df_total = pd.DataFrame()
pre_df_total = pd.DataFrame()
recall_df_total = pd.DataFrame()
session_name = "train"
type_data = "HFI_NI_PROM_NEW"
folder_number = 0
for strain_name in strain_names:
    accuracy_file = "/home/mddo/stage/M2S4/output/{}/accuracy/{}/accuracy_{}_{}.csv".format(strain_name, session_name, type_data,folder_number)
    # create_folder("/home/mddo/stage/M2S4/images/{}".format(strain_name))
    # create_folder("/home/mddo/stage/M2S4/images/{}/acc_precision_recall/".format(strain_name))
    # save_figure = "/home/mddo/stage/M2S4/images/{}/acc_precision_recall/{}_{}.png".format(strain_name,type_data, folder_number)
    accuracy_df = pd.read_csv(accuracy_file)
    accuracy_df.columns = ["forest","accuracy","precision","recall","fscrore","total_tree"]

    accuracy = accuracy_df["accuracy"]
    precision = accuracy_df["precision"]
    total_tree = accuracy_df["total_tree"]
    acc_df_total["acc_{}".format(strain_name)] = accuracy_df["accuracy"]
    acc_df_total["prec_{}".format(strain_name)] = accuracy_df["precision"]
    acc_df_total["recall_{}".format(strain_name)] = accuracy_df["recall"]
# array = ["accuracy","precision","recall"]
# for key in array:
# if key == "accuracy":
#     df_plot = acc_df_total
# elif key == "precision":
#     df_plot = pre_df_total
# else:
#     df_plot = recall_df_total
ax = plt.gca()
sns.set_style("whitegrid")

acc_df_total.plot(kind='line',ax=ax)
# accuracy_df.plot(kind='line',y='precision', color='red', ax=ax)
# accuracy_df = accuracy_df.drop(columns = ["total_tree"])
# accuracy_df.plot(kind="line")

fig = matplotlib.pyplot.gcf()
fig.set_size_inches(18.5, 8.5)
plt.rcParams["figure.figsize"] = (10,2)
# plt.title("{} prediction accuracy and precision".format(strain_name), size = 20)
plt.xticks(size = 14)
plt.yticks(size = 14)
plt.xlabel("Iteration", size = 18)
plt.ylabel("", size = 18)
# plt.xlim(0,60)

# plt.savefig("/home/mddo/stage/M2S4/images/{}_{}_{}_{}.png".format(key, session_name,type_data, folder_number))
# plt.clf()

# %%
df = pd.read_csv("/home/mddo/stage/M2S4/output/FY/diploid_/diploid_0/df/normal_raw.csv")

# %%
sns.boxplot(x = "insertion_index", y = "label", data = df)

# %%
#####GROUP BARPLOT FALSE POSITVIVE AND FALSE NEGATIVE
session_name = "train"
type_df = "normal_KNN"
labels = ["FY"]
# labels = ["FY","ABP","ACF","ACN","ACP","ADD","AND","APH","AVI","BBQ","BHH","BMK","CCD","CGQ","CHM","CIB","CLG","CNM","CNT","CPG","Sigma"]
FP_means = []
FN_means = []
for strain in labels:
    FP_file = "/home/mddo/stage/M2S4/output/{}/error/{}/{}_0_FP.csv".format(strain, session_name, type_df)
    FN_file = "/home/mddo/stage/M2S4/output/{}/error/{}/{}_0_FN.csv".format(strain, session_name, type_df)

    fp_df = pd.read_csv(FP_file)
    fp_df.columns = ["orf","freq"]
    count_fp,_ = fp_df[fp_df["freq"] >= 10].shape
    FP_means.append(count_fp)

    fn_df = pd.read_csv(FN_file)
    fn_df.columns = ["orf","freq"]
    count_fn,_ = fn_df[fn_df["freq"] >= 10].shape
    FN_means.append(count_fn)

print(FP_means)
print(FN_means)
x = np.arange(len(labels))  # the label locations
width = 0.35  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x - width/2, FP_means, width, label='FP')
rects2 = ax.bar(x + width/2, FN_means, width, label='FN')

# Add some text for labels, title and custom x-axis tick labels, etc.
# ax.set_title('Number of FP and FN for each strain')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()


def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')


autolabel(rects1)
autolabel(rects2)

fig.set_size_inches(18.5, 8.5)
plt.rcParams["figure.figsize"] = (10,2)
# plt.xticks(rotation=45, ha='right')
plt.xticks(size = 14)
plt.yticks(size = 14)
# plt.xlim([-1, X.shape[1]])
plt.xlabel("Strains", size = 20)
plt.ylabel('Number of gene', size = 20)
plt.show()

# %%
