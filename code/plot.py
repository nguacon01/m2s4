#%%
import numpy as numpy
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import confusion_matrix
from helper_functions import train_test_split
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
report_normal_path = "/home/mddo/stage/M2S4/output/accuracy/test/accuracy_normal.csv"

# %%
df_acc_HFI_NI_PROM = pd.read_csv(report_HFI_NI_PROM_path).sort_values(by=["accuracy"])
df_acc_HFI_NI_PROM.reset_index(inplace=True)
print(df_acc_HFI_NI_PROM)

# %%
df_acc_normal_test = pd.read_csv(report_normal_path).sort_values(by=["accuracy"])
df_acc_normal_test.reset_index(inplace=True)

# %%
merge_df_test = pd.concat([df_acc_normal["accuracy"],df_acc_HFI_NI_PROM["accuracy"]], axis = 1, join="inner")
merge_df_test.columns = ["acc_normal","acc_HFI_NI_PROM"]
#%%
print(merge_df_test)

# %%
sns.set(style="whitegrid")
df_acc_normal_test["accuracy"].plot(kind="line")
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(18.5, 10.5)
plt.rcParams["figure.figsize"] = (20,2)

# %%
df_full = pd.read_csv("/home/mddo/stage/M2S4/output/FY/haploid/df/train/normal.csv")



# %%
