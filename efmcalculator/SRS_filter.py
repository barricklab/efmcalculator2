import pandas as pd

def filter_redundant(df):
<<<<<<< HEAD
    df.sort_values(by=['length', 'sequence'], ascending=[False, True], inplace=True)
    df = df.reset_index(drop=True)
    i = 0
    while i < len(df):
        deleted = False
        # creates new df with only repeats that contain the sequence of the row
        similar = df[df['sequence'].str.contains(df.loc[i, 'sequence'])]
        # if similar contains a repeat other than the repeat we're looking at
        if len(similar) > 1:
            # resets the index values of the similar df
            similar = similar.reset_index(drop=True)

            x = 0
            while x < len(similar) and not deleted:
                # check if this is the repeat we were looking at
                if similar.loc[x, 'length'] > df.loc[i, 'length']:
                    # check if they have the same number of positions
                    if similar.loc[x, 'occurrences'] == df.loc[i, 'occurrences']:
                        # delete the row we were looking at
                        df.drop([i], inplace=True)
                        df = df.reset_index(drop=True)
                        i -= 1
                        deleted = True
                        print("redundant")
                x += 1
=======
    i = 0
    length = len(df)
    while i < length:
        # creates new df with only repeats that contain the sequence of the row
        similar = df[df['Sequence'].str.contains(df.loc[i, 'Sequence'])]
        print(df.loc[i, 'Sequence'])
        # resets the index values of the similar df
        similar = similar.reset_index(drop=True)
        x = 0
        count = df.loc[i, 'Occurrences']
        y = 0
        deleted = False
        while y < len(similar) and not deleted:
            if similar.loc[y, 'Sequence'] == df.loc[i, 'Sequence']:
                if count < len(similar) and count == similar.loc[x + count, 'Occurrences']:
                    while x < count:
                        print("drop", df.loc[i, 'Sequence'])
                        df.drop([i], inplace=True)
                        x += 1
                        i += 1
                        deleted = True
            y += 1

>>>>>>> c2caf551d6125e7854ea443fccf8bc983bad780a
        i += 1
    return df

    

<<<<<<< HEAD
def consolidate_repeats(df):
    new_df = None
    grouped = df.groupby('sequence')
    print(grouped.head())
    return new_df
=======
start = time.time()

file = r"G:\Other computers\My Computer\School\EvoStab\Computation\psl9Sreya_3.csv"

df = pd.read_csv(file)


# Add occurrences column to df
occurrences = df['Sequence'].value_counts()
df['Occurrences'] = df['Sequence'].map(occurrences)
df.sort_values(by=['Size', 'Sequence'], ascending=[True, True], inplace=True)

df = filter_redundant(df)
print(df)
df.to_csv(r"G:\Other computers\My Computer\School\EvoStab\Computation\psl9Sreya_3_filter_5.csv", index=False)

end = time.time()
print("There are " + str(len(df)) + " repeated sequences")
print(end-start)
>>>>>>> c2caf551d6125e7854ea443fccf8bc983bad780a
