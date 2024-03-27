def filter_redundant(df):
    df.sort_values(by=["length", "sequence"], ascending=[False, True], inplace=True)
    df = df.reset_index(drop=True)
    i = 0
    while i < len(df):
        deleted = False
        # creates new df with only repeats that contain the sequence of the row
        similar = df[df["sequence"].str.contains(df.loc[i, "sequence"])]
        # if similar contains a repeat other than the repeat we're looking at
        if len(similar) > 1:
            # resets the index values of the similar df
            similar = similar.reset_index(drop=True)

            x = 0
            while x < len(similar) and not deleted:
                # check if this is the repeat we were looking at
                if similar.loc[x, "length"] > df.loc[i, "length"]:
                    # check if they have the same number of positions
                    if similar.loc[x, "occurrences"] == df.loc[i, "occurrences"]:
                        # delete the row we were looking at
                        df.drop([i], inplace=True)
                        df = df.reset_index(drop=True)
                        i -= 1
                        deleted = True
                        print("redundant")
                x += 1
        i += 1
    return df
