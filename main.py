from help_function import *
def main():
    df = load_df("iris.csv")
    train_df,test_df = train_test_split(df, 0.2)
    print(test_df.head())
if __name__ == "__main__":
    main()