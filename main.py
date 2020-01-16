from help_function import *
def main():
    df = load_df("iris.csv")
    train_df, test_df = train_test_split(df, 0.2)
    data = train_df.values
    classi = classify_data(data)
    print(classi)
if __name__ == "__main__":
    main()