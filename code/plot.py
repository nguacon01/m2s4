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
