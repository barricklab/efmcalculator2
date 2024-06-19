def filter_redundant(df):
    df.sort_values(by=["length", "sequence"], ascending=[False, True], inplace=True)
    df = df.reset_index(drop=True)
    i = 0
    while i < len(df):
        deleted = False

        # Delete rows with overlapping repeats and only 2 occurrences (not a real repeat)
        if df.loc[i, "occurrences"] == 2:
            first_pos = df.loc[i, "positions"][0]
            second_pos = df.loc[i, "positions"][1]
            length = df.loc[i, "length"]
            if (
                abs(second_pos[0] - first_pos[0]) < length
                or abs(second_pos[1] - first_pos[1]) < length
            ):
                df.drop([i], inplace=True)
                df = df.reset_index(drop=True)
                i -= 1
                deleted = True

        # Corrects output of SSR with length >= 6 and <= 15

        if not deleted:
            if df.loc[i, "occurrences"] > 1:
                if df.loc[i, "positions"][1][2] == "SSR":
                    # print("detected")
                    repeat_count = int(df.loc[i, "positions"][1][3])
                    start_pos = df.loc[i, "positions"][1][0]
                    end_pos = df.loc[i, "positions"][1][1]
                    repeat_length = df.loc[i, "length"]
                    mut_rate = df.loc[i, "positions"][1][4]
                    new_positions = [
                        (start_pos, end_pos, "SSR", repeat_count, mut_rate)
                    ]
                    # replaces positions column in df with correct information
                    df.at[i, "positions"] = new_positions
                    df.at[i, "occurrences"] = 1

            # Delete redundant SRS and RMD repeats, and Keep SSR with lowest length

            # creates new df with only repeats that contain the sequence of the row
            similar = df[
                df["sequence"].str.contains(df.loc[i, "sequence"])
                & df["positions"].apply(
                    lambda x: x[0][2] == df.loc[i, "positions"][0][2]
                )
            ]

            # if similar contains a repeat other than the repeat we're looking at
            if len(similar) > 1:
                # resets the index values of the similar df
                similar = similar.reset_index(drop=True)

                # Delete redundant SRS repeats

                x = 0
                while x < len(similar) and not deleted:
                    # Delete redundant SRS repeats
                    if df.loc[i, "positions"][0][2] in {"SRS", "RMD"}:
                        # check if longer repeat is found
                        if similar.loc[x, "length"] > df.loc[i, "length"]:
                            # check if they have the same number of positions
                            if (
                                similar.loc[x, "occurrences"]
                                == df.loc[i, "occurrences"]
                            ):
                                # delete the row we were looking at
                                df.drop([i], inplace=True)
                                df = df.reset_index(drop=True)
                                i -= 1
                                deleted = True

                    # Keep SSR with lowest length

                    elif df.loc[i, "positions"][0][2] == "SSR":
                        # make sure both SSRs are starting from the same position
                        if df.loc[i, "positions"][0][0] == df.loc[x, "positions"][0][0]:
                            # Compare lengths of both SSRs
                            if df.loc[i, "length"] > df.loc[x, "length"]:
                                df.drop([i], inplace=True)
                                df = df.reset_index(drop=True)
                                i -= 1
                                deleted = True
                    x += 1
        i += 1
    return df
