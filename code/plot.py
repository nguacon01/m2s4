#%%
import numpy as numpy
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import confusion_matrix
from helper_functions import train_test_split

df = pd.read_csv("output/FY/dataframe.csv")
df.shape
train_df, test_df = train_test_split(df, 0.2)
data = test_df['label'].value_counts()
print(data)
data.plot(kind="bar",color = ["blue","orange"])
plt.xticks(rotation=0,size=16)
plt.yticks(size=14)
plt.ylabel("numbers of gene",size=14)
plt.show()
ax = sns.barplot(df['label'].value_counts())



# ax = sns.scatterplot(x="insertion_index", y="NI_index", hue="label", data=df)

# sns.heatmap([train_df.hits_count,train_df.read_count],annot=True,fmt="d")

# %%
sns.scatterplot(y='HFI_normalized', x='hits_count', data = df, hue='label')

# %%
df = pd.read_csv("output/output_predictions_epoche_0.csv")
confusion_matrix = pd.crosstab(df['label'],df['predictions'], rownames = ['Actual'], colnames=['Predict'])
sns.heatmap(confusion_matrix,annot=True,annot_kws={"size": 16,"color":"white"})
plt.show()



# %%
df['hits_count'].plot(kind='bar')

# %%
ver2_df = pd.read_csv("/home/mddo/stage/M2S4/output/FY/dataframe.csv")

print(ver2_df.head(10))# %%


# %%


# %%
sns.scatterplot(y='insertion_index', x='ratio_hits_prom', data = ver2_df, hue='label')
plt.yticks(size = 14)
plt.xticks(size = 14)
plt.xlabel("ratio hits prom", size = 16)
plt.ylabel("hits count in 100 bp prom", size = 16)
plt.rcParams["figure.figsize"] = (30,4)

# %%
dataFrame = pd.read_csv("/home/mddo/stage/M2S4/output/FY/accuracy_para.csv")
# %%
data = pd.DataFrame(dataFrame["accuracy"])
# %%
sns.lineplot(data=data, palette="tab10", linewidth=2.5)
plt.rcParams["figure.figsize"] = (30,4)
sns.set(style="whitegrid")

# %%
