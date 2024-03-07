import time
import pandas as pd


def filter_redundant(df):
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

        i += 1
    return df


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
