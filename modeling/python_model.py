from statsmodels.gam.api import GLMGam, CyclicCubicSplines
import csv
import polars

def main():
    dataframe = polars.read_csv(source="input.csv")
    data_of_interest = dataframe.select(["Mutation_Rate", "RBP_Length", "TBD_length"])
    spline = CyclicCubicSplines(x=data_of_interest, df=[3,3, 3])

    gam = GLMGam.from_formula('Mutation_Rate ~ RBP_Length + TBD_length', data=data_of_interest,
                                  smoother=spline)

    res = gam.fit()

    predict_data = polars.dataframe({"RBP_Length": 0,
                                     "TBD_length": 0})

    print(res.summary())


if __name__ == "__main__":
    main()