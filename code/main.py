#%%
from help_function import *
from helper_functions import *
from random_forest import *
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn import metrics
import glob
import json
def main():

    #CREATE ORIGINAL DATA
    # strain_names = ["Sigma","CNT","CPG"]
    strain_name = "Sigma"
    # # #---------------merge data file--------------#
    merge_df(strain_name)

    ##--------------------# BEGIN TESTING_SESSION#--------------------##
    # i = 0
    # type_data = "HFI_NI_PROM_KNN_{}".format(i)
    # strain_name = "CNT"
    
    # test_df_path = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/df/test/HFI_NI_PROM_KNN.csv".format(strain_name,i)
    # testing_RF(test_df_path, type_data, strain_name)

    ##--------------------# END TESTING_SESSION#--------------------##


    


if __name__ == "__main__":
    main()




# %%
