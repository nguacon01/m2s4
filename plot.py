#%%
import numpy as numpy
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv("output/df_df.csv")

print(df['label'].value_counts())
data = df['label'].value_counts()
print(data)
data.plot(kind="bar")
plt.ylabel("numbers of gene")
plt.show()
# ax = sns.barplot(df['label'].value_counts())



# ax = sns.scatterplot(x="insertion_index", y="NI_index", hue="label", data=df)

# sns.heatmap([train_df.hits_count,train_df.read_count],annot=True,fmt="d")

# %%
sns.scatterplot(x='nonCoding_windows', y='insertion_index', data = df, hue='label')

# %%
