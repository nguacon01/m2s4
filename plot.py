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
sns.scatterplot(y='hits_count_pro', x='hits_count', data = df, hue='label')

# %%
df = pd.read_csv("output/output_predictions_epoche_0.csv")
confusion_matrix = pd.crosstab(df['label'],df['predictions'], rownames = ['Actual'], colnames=['Predict'])
sns.heatmap(confusion_matrix,annot=True,annot_kws={"size": 16,"color":"white"})
plt.show()



# %%
df['hits_count'].plot(kind='bar')

# %%
