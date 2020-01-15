import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import random

def load_df(file_name):
    df = pd.read_csv(file_name)
    return df

def train_test_split(df, test_size):
    if isinstance(test_size, float):
        test_size = round(test_size * len(df))
    indices = df.index.tolist()
    #choose random indices
    test_indices = random.sample(population=indices,k=test_size)
    test_df = df.loc[test_indices]
    train_df = df.drop(test_indices)
    return train_df, test_df