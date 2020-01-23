#%%
import numpy as np
import random
import pandas as pd
from helper_functions import train_test_split

df = pd.read_csv("iris.csv")
train_df, test_df = train_test_split(df, 0.2)

#randomly create sub dataframes based on train_df
#input: train_df
#n_bootstrap: number of rows within sub dataframe
def bootstrapping(train_df, n_bootstrap):
    bootstrap_indices = np.random.randint(low=0,high=len(train_df), size=n_bootstrap)
    df_bootstrapped = train_df.iloc[bootstrap_indices]

    return df_bootstrapped

# %%
bs_df = bootstrapping(train_df, n_bootstrap = 100).duplicated()

# %%


# %%
