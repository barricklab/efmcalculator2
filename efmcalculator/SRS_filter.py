import time
import pandas as pd


def filter_redundant(df):
    memory = [""]
    i = 0
    while i < len(df):
        # creates new df with only repeats that contain the sequence of the row
        similar = df[df['Sequence'].str.contains(df.loc[i, 'Sequence'])]
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
                        print("drop", df.loc[i+x, 'Sequence'])
                        df.drop([i+x], inplace=True)
                        x += 1
                        deleted = True
                    i -= count
            y += 1
            df = df.reset_index(drop=True)
        i += count
    return df


start = time.time()

file = "G:\Other computers\My Computer\School\EvoStab\Computation\psl9Sreya_3.csv"

df = pd.read_csv(file)


# Add occurrences column to df
occurrences = df['Sequence'].value_counts()
df['Occurrences'] = df['Sequence'].map(occurrences)
df.sort_values(by=['Size', 'Sequence'], ascending=[True, True], inplace=True)

df = filter_redundant(df)
print(df)
df.to_csv("G:\Other computers\My Computer\School\EvoStab\Computation\psl9Sreya_3_filtered_3.csv", index=False)

end = time.time()
print("There are " + str(len(df)) + " repeated sequences")
print(end-start)
