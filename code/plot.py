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
# %%
sns.lineplot(data=pd.DataFrame(dataFrame["accuracy_HFI_NI"]), palette="tab10", linewidth=2.5, color = "orange",markers="o")

# plt.rcParams["figure.figsize"] = (30,4)
# sns.set(style="whitegrid")
# dataFrame.plot(kind="line")
plt.rcParams["figure.figsize"] = (30,4)
sns.set(style="whitegrid")

# %%
#PLOT TRAINING ACCURACY
df_HFI_NI = pd.read_csv("/home/mddo/stage/M2S4/output/accuracy/train/accuracy_HFI_NI.csv")
df_HFI_NI_PROM = pd.read_csv("/home/mddo/stage/M2S4/output/accuracy/train/accuracy_HFI_NI_PROM.csv")
df_normal = pd.read_csv("/home/mddo/stage/M2S4/output/accuracy/train/accuracy_normal.csv")
df_HFI_PROM = pd.read_csv("/home/mddo/stage/M2S4/output/accuracy/train/accuracy_HFI_PROM.csv")
df_NI_PROM = pd.read_csv("/home/mddo/stage/M2S4/output/accuracy/train/accuracy_NI_PROM.csv")

#%%
print(df_HFI_NI["accuracy"])
# %%

df_HFI_NI_acc = pd.DataFrame(df_HFI_NI["accuracy"]).sort_values(by=["accuracy"])
df_HFI_NI_acc.reset_index(inplace=True)
df_HFI_NI_PROM_acc = pd.DataFrame(df_HFI_NI_PROM["accuracy"]).sort_values(by=["accuracy"])
df_HFI_NI_PROM_acc.reset_index(inplace=True)
df_normal_acc = pd.DataFrame(df_normal["accuracy"]).sort_values(by=["accuracy"])
df_normal_acc.reset_index(inplace=True)
df_HFI_PROM_acc = pd.DataFrame(df_HFI_PROM["accuracy"]).sort_values(by=["accuracy"])
df_HFI_PROM_acc.reset_index(inplace=True)
df_NI_PROM_acc = pd.DataFrame(df_NI_PROM["accuracy"]).sort_values(by=["accuracy"])
df_NI_PROM_acc.reset_index(inplace=True)
#%%
merge_df = pd.concat([df_normal_acc["accuracy"], df_HFI_NI_PROM_acc["accuracy"]],axis=1, join='inner')
merge_df.columns = ["df_normal_acc", "df_HFI_NI_PROM_acc"]
# %%
sns.set(style="whitegrid")
merge_df.plot(kind="line")
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(18.5, 10.5)
plt.rcParams["figure.figsize"] = (20,2)
# %%

#PLOT TESTING ACCURACY
report_HFI_NI_PROM_path = "/home/mddo/stage/M2S4/output/accuracy/test/accuracy_HFI_NI_PROM.csv"
report_HFI_NI_path = "/home/mddo/stage/M2S4/output/accuracy/test/accuracy_HFI_NI.csv"
report_HFI_PROM_path = "/home/mddo/stage/M2S4/output/accuracy/test/accuracy_HFI_PROM.csv"
report_NI_PROM_path = "/home/mddo/stage/M2S4/output/accuracy/test/accuracy_NI_PROM.csv"
report_normal_path = "/home/mddo/stage/M2S4/output/accuracy/test/accuracy_normal.csv"

# %%
df_acc_HFI_NI_PROM = pd.read_csv(report_HFI_NI_PROM_path).sort_values(by=["accuracy"])
df_acc_HFI_NI_PROM.reset_index(inplace=True)
# %%
df_acc_HFI_NI = pd.read_csv(report_HFI_NI_path).sort_values(by=["accuracy"])
df_acc_HFI_NI.reset_index(inplace=True)
# %%
df_acc_HFI_PROM = pd.read_csv(report_HFI_PROM_path).sort_values(by=["accuracy"])
df_acc_HFI_PROM.reset_index(inplace=True)
# %%
df_acc_NI_PROM = pd.read_csv(report_NI_PROM_path).sort_values(by=["accuracy"])
df_acc_NI_PROM.reset_index(inplace=True)
# %%
df_acc_normal = pd.read_csv(report_normal_path).sort_values(by=["accuracy"])
df_acc_normal.reset_index(inplace=True)

# %%
merge_df_test = pd.concat([df_acc_HFI_NI_PROM["accuracy"],df_acc_HFI_NI["accuracy"],df_acc_HFI_PROM["accuracy"],df_acc_NI_PROM["accuracy"]], axis = 1, join="inner")
merge_df_test.columns = ["acc_HFI_NI_PROM","acc_HFI_NI","acc_HFI_PROM","acc_NI_PROM"]
#%%
print(merge_df_test)

# %%
sns.set(style="whitegrid")
merge_df_test.plot(kind="line")
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(18.5, 10.5)
plt.rcParams["figure.figsize"] = (20,2)

# %%
df_full = pd.read_csv("/home/mddo/stage/M2S4/output/FY/haploid/df/train/normal.csv")



# %%
seed(1)
np.random.randint(low=0,high=2494, size=1800)

# %%


# %%
seed(1)
np.random.randint(low=0,high=2494, size=1800)

# %%


# %%
df = pd.read_csv("/home/mddo/stage/M2S4/output/FY/diploid_/diploid_0/df/HFI_NI_PROM.csv")

# %%
print(df)

# %%
NI = pd.read_csv("/home/mddo/stage/M2S4/output/FY/diploid_/diploid_0/diplo_NI_ratio_haplo_diplo.out", sep=" ")
print(NI)

# %%

"""
REFRORM STRUCTURE OF ACCURACY FILE
NAME OF FOREST _ ACCURACY _ TOTAL NUMBER OF TREE IN THE FOREST

"""
acc_file_path = "/home/mddo/stage/M2S4/output/accuracy/test/accuracy_HFI_NI_PROM_nan.csv"
df = pd.read_csv(acc_file_path)
print(df)

# %%
forest_name = df["forest"]
print(forest_name)

# %%
total_num_tree_arr = []
forest_name_arr = []
for fn in forest_name:
    forest_name = "{}.0".format(fn)
    forest_name_arr.append(forest_name)
    forest_path = "/home/mddo/stage/M2S4/output/forest/HFI_NI_PROM_nan/{}.json".format(forest_name)
    with open(forest_path) as content:
        json_data = json.load(content)
        total_numner_tree = len(json_data)
    total_num_tree_arr.append(total_numner_tree)
        

# %%
print(df.shape)
print(len(total_num_tree_arr))

# %%
df["total_tree"] = total_num_tree_arr
df["forest"] = forest_name_arr

# %%
print(df)
print(df.columns)
df.drop(df.columns.difference(["forest","accuracy","total_tree"]), 1, inplace = True)


# %%
df.to_csv(acc_file_path, index = False)

"""
FINISH REFORM ACCURACY FILE
"""
# %%

"""
CREATE FIGURE OF TRAININGS WITH 2 DATAFRAME: WITH NaN VALUE AND non NaN VALUE
"""

df_acc_train_HFI_NI_PROM_path = "/home/mddo/stage/M2S4/output/accuracy/train/accuracy_HFI_NI_PROM.csv"
df_acc_train_HFI_NI_PROM_nan_path = "/home/mddo/stage/M2S4/output/accuracy/train/accuracy_HFI_NI_PROM_nan.csv"
df_acc_test_HFI_NI_PROM_path = "/home/mddo/stage/M2S4/output/accuracy/test/accuracy_HFI_NI_PROM.csv"
df_acc_test_HFI_NI_PROM_nan_path = "/home/mddo/stage/M2S4/output/accuracy/test/accuracy_HFI_NI_PROM_nan.csv"

df_acc_train_HFI_NI_PROM = pd.read_csv(df_acc_train_HFI_NI_PROM_path)
df_acc_train_HFI_NI_PROM.sort_values(by=["accuracy"])
df_acc_train_HFI_NI_PROM.reset_index(inplace=True)
df_acc_train_HFI_NI_PROM_nan = pd.read_csv(df_acc_train_HFI_NI_PROM_nan_path)
df_acc_train_HFI_NI_PROM_nan.sort_values(by=["accuracy"])
df_acc_train_HFI_NI_PROM_nan.reset_index(inplace=True)
df_acc_test_HFI_NI_PROM = pd.read_csv(df_acc_test_HFI_NI_PROM_path)
df_acc_test_HFI_NI_PROM.sort_values(by=["accuracy"])
df_acc_test_HFI_NI_PROM.reset_index(inplace=True)
df_acc_test_HFI_NI_PROM_nan = pd.read_csv(df_acc_test_HFI_NI_PROM_nan_path)
df_acc_test_HFI_NI_PROM_nan.sort_values(by=["accuracy"])
df_acc_test_HFI_NI_PROM_nan.reset_index(inplace=True)

#%%
merge_df = pd.concat([df_acc_train_HFI_NI_PROM_nan["accuracy"],df_acc_test_HFI_NI_PROM_nan["accuracy"]], axis = 1, join="inner")
merge_df.columns = ["acc_train_HFI_NI_PROM_nan","acc_test_HFI_NI_PROM_nan"]
#%%
sns.set(style="whitegrid")
merge_df.plot()
# sns.lineplot(x = "total_tree", y = "accuracy", data = merge_df)
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(18.5, 10.5)
plt.rcParams["figure.figsize"] = (20,2)


# %%
"""
END
"""
#%%

predictions_HFI_NI_PROM_train = pd.read_csv("/home/mddo/stage/M2S4/output/predictions/train/HFI_NI_PROM/predictions_18_4_10_2106_92.0.csv")
label_actu = predictions_HFI_NI_PROM_train["label"]
label_pred = predictions_HFI_NI_PROM_train["predictions"]
confusion_matrix = confusion_matrix(label_actu,label_pred)
sns.heatmap(confusion_matrix)

# %%
df = pd.read_csv("/home/mddo/stage/M2S4/data/FY/final_annot.csv", sep = ";")

# %%
print(df["systematic_name"])

# %%
df = df.drop(columns = ["AnnotSGD", "AnnotDowell", "NoteFinal"])

# %%
print(df)

# %%
drop_index = df[df["AnnotFinal"] == "decrease"].index
print(drop_index)

#%%
df.drop(drop_index, inplace = True)

# %%
drop_index = df[df["AnnotFinal"] == "decrease_galactose"].index
print(drop_index)
df.drop(drop_index, inplace = True)
df.shape

# %%
drop_index = df[df["AnnotFinal"] == "Divergent"].index
print(drop_index)
df.drop(drop_index, inplace = True)
df.shape

# %%
df.to_csv("/home/mddo/stage/M2S4/data/FY/final_annot.csv", index = False)

# %%
df = df.replace("No","non_ess")

# %%
df = df.replace("Yes","ess")

# %%
df.to_csv("/home/mddo/stage/M2S4/data/FY/final_annot.csv", index = False)

# %%
df = pd.read_csv("/home/mddo/stage/M2S4/output/FY/diploid_/diploid_0/df/HFI_NI_PROM.csv")

# %%
df.info()

# %%


# %%
df = pd.read_csv("/home/mddo/stage/M2S4/output/predictions/train/HFI_NI_PROM_zerofill_0/predictions_forest_19_5_10_1789_92.0.csv")
plt.figure(figsize=(10,10))
confusion_matrix = pd.crosstab(df['label'],df['predictions'], rownames = ['Actual'], colnames=['Predict'])
ax = sns.heatmap(confusion_matrix,
            annot=True,
            annot_kws={"size": 22,},
            fmt='g',
            vmin=0, vmax=600,
            linewidths=.5,
            cbar=False)
plt.xticks(size=20)
plt.yticks(size=20)
plt.xlabel("Predict",size=14)
plt.ylabel("Actual", size=14)

bottom, top = ax.get_ylim()
ax.set_ylim(bottom + 0.5, top - 0.5)

# plt.rcParams.update({'font.size': 14})
# plt.show()# %%
plt.savefig("/home/mddo/stage/M2S4/images/confusion_matrix_train_zerofill.png")


# %%
y_true = df["label"]
y_pred = df["predictions"]
precision_recall_fscore_support(y_true, y_pred, average='macro')

# %%
df = pd.read_csv("/home/mddo/stage/M2S4/output/predictions/train/HFI_NI_PROM_zerofill_0/predictions_forest_19_5_10_1789_92.0.csv")

# %%
y_true = df["label"]
y_pred = df["predictions"]
metric = precision_recall_fscore_support(y_true, y_pred, )


# %%
print(metric)

# %%
