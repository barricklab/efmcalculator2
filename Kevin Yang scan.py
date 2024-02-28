import regex as re
import pandas as pd
from Bio import SeqIO


# Finds the positions of a given repeat
def positions(position, middle):
    if not middle:
        for z in re.finditer(n, my_seq, overlapped=True):
            position.append(z.start() + 1)
    else:
        for z in re.finditer(n, my_seq2, overlapped=True):
            position.append(z.start() + 1)


# Finds the distances between repeats
def distances(position, distance):
    z = 1
    circ_dist = []
    while z < len(position):
        d = position[z] - position[z - 1]
        distance.append(d)
        d2 = (position[z - 1] + len(my_seq)) - (position[z])
        circ_dist.append(d2)
        z += 1
    d = 0
    while d < len(distance):
        if distance[d] > circ_dist[d]:
            distance[d] = circ_dist[d]
        d += 1


# Adds repeats to the dataframe
def add(table, position, distance, number):
    temp_df = pd.DataFrame({'Sequence': [str(n)],
                            'Length': [length],
                            'Occurrences': [number],
                            'Positions': [position],
                            'Distances': [distance],
                            })
    table = table._append(temp_df, ignore_index = True)
    return table


def unique(table):
    u = 0
    while u < len(table):
        if n in table.loc[u]['Sequence']:
            return False
        u += 1
    return True


def update(table, pos, index, locations):
    # updates positions
    locations.append(pos + 1)
    # updates distances
    dis = table.loc[index[0], 'Distances']
    dis.append(0)
    # gets correct distances
    dist = []
    distances(locations, dist)
    j = 0
    # puts correct distances in table
    while j < len(dis):
        dis[j] = dist[j]
        j += 1
    # finds previous number of occurrences
    amount = table.loc[index[0], 'Occurrences']
    # updates number of occurrences
    table.loc[index[0], 'Occurrences'] = amount + 1


# Enter sequence using fasta file
"""
filename = "psl9.fasta"
sequence = SeqIO.read(filename, 'fasta')
my_seq = str(sequence.seq).upper()
"""

# Enter sequence as a string
sequence = input("what is the sequence?")
my_seq = sequence.upper()

data = {'Sequence': [[]],
        'Length': [[]],
        'Occurrences': [[]],
        'Positions': [[]],
        'Distances': [[]],
        }
df = pd.DataFrame(data)


# Finds all repeats starting from beginning of sequence
i = 0
length = 6

# Adds all 6 bp repeats to dataframe
print("finding 6 bp repeats")
while i + length - 1 < len(my_seq):
    pos = []
    dist = []
    n = my_seq[i:i + length]
    count = len(re.findall(n, my_seq, overlapped=True))

    if unique(df) and count >= 2:
        # find positions of all repeats
        positions(pos, False)

        # find distances between repeats
        distances(pos, dist)

        # adds repeat to dataframe
        df = add(df, pos, dist, len(pos))
    print("i = " + str(i))
    i += 1

# finds larger repeats
print("finding larger repeats")
i = 1
while i < len(df):
    length = df.loc[i, 'Length'] + 1
    j = 0
    # j and k used to do pairwise comparison of all positions
    while j < df.loc[i, 'Occurrences'] - 1:
        k = j + 1
        while k < df.loc[i, 'Occurrences']:
            x = df.loc[i, 'Positions'][j] - 1
            y = df.loc[i, 'Positions'][k] - 1
            # if a larger repeat is found, and they are not overlapping
            if my_seq[x:x+length] == my_seq[y:y + length] and x+length <= y:
                n = my_seq[x:x+length]
                # If a new row with the larger repeat has already been created
                if len(df[df['Sequence'] == my_seq[x:x+length]]) > 0:
                    # finds index position of repeat in df
                    index = df.loc[df['Sequence'] == my_seq[x:x + length]].index
                    # finds previous positions and stores it in new list
                    locations = df.loc[index[0], 'Positions']
                    # if y position is not included in row
                    if locations.count(y + 1) == 0:
                       update(df, y, index, locations)
                    # if x position is not included in row
                    if locations.count(x + 1) == 0:
                        update(df, x, index, locations)
                # if no new row was created, create a new row
                else:
                    dist = []
                    distances([x+1, y+1], dist)
                    df = add(df, [x + 1, y + 1], dist, 2)
            k += 1
        j += 1
    i += 1

# Finds all repeats starting from the middle of the sequence
my_seq2 = my_seq[int(len(my_seq)/2):] + my_seq[:int(len(my_seq)/2)]
df2 = pd.DataFrame(data)

i = 0
length = 6

# Adds all 6 bp repeats from middle of the sequence to dataframe
print("finding 6 bp repeats from the middle")
while i + length - 1 < len(my_seq2):
    pos = []
    dist = []
    n = my_seq2[i:i + length]
    count = len(re.findall(n, my_seq2, overlapped=True))

    if unique(df2) and count >= 2:
        # find positions of all repeats
        positions(pos, True)
        pos.sort()

        # find distances between repeats
        distances(pos, dist)

        # adds repeat to dataframe
        df2 = add(df2, pos, dist, len(pos))
    print("i2 = " + str(i))
    i += 1

# finds larger repeats from middle of the sequence
print("finding larger repeats from middle")
i = 1

while i < len(df2):
    length = df2.loc[i, 'Length'] + 1
    j = 0
    # j and k used to do pairwise comparison of all positions
    while j < df2.loc[i, 'Occurrences'] - 1:
        k = j + 1
        while k < df2.loc[i, 'Occurrences']:
            x = df2.loc[i, 'Positions'][j] - 1
            y = df2.loc[i, 'Positions'][k] - 1
            # if a larger repeat is found, and they are not overlapping
            if my_seq2[x:x+length] == my_seq2[y:y + length] and x+length <= y:
                n = my_seq2[x:x+length]
                # If a new row with the larger repeat has already been created
                if len(df2[df2['Sequence'] == my_seq2[x:x+length]]) > 0:
                    # finds index position of repeat in df
                    index = df2.loc[df2['Sequence'] == my_seq2[x:x+length]].index
                    # finds previous positions and stores it in new list
                    locations = df2.loc[index[0], 'Positions']
                    # check if y position is already included in the row
                    if locations.count(y + 1) == 0:
                        update(df2, y, index, locations)
                    # check if x position is already included in the row
                    if locations.count(x + 1) == 0:
                        update(df2, x, index, locations)
                # if no new row was created, create a new row
                else:
                    dist = []
                    distances([x+1, y+1], dist)
                    df2 = add(df2, [x + 1, y + 1], dist, 2)
            k += 1
        j += 1
    i += 1


# deletes first row (which is empty) in both tables
df.drop([0], inplace=True)
df = df.reset_index(drop=True)
df2.drop([0], inplace=True)
df2 = df2.reset_index(drop=True)


# changes positions in df2 to be their positions when starting from the beginning
i = 0
while i < len(df2):
    j = 0
    while j < df2.loc[i, 'Occurrences']:
        if len(my_seq) % 2 == 0:
            if df2.loc[i, 'Positions'][j] <= int(len(my_seq2) / 2):
                df2.loc[i, 'Positions'][j] += int(len(my_seq2) / 2)
            else:
                p = df2.loc[i, 'Positions'][j]
                df2.loc[i, 'Positions'][j] = p - int(len(my_seq2) / 2)
        else:
            if df2.loc[i, 'Positions'][j] <= int(len(my_seq2) / 2):
                df2.loc[i, 'Positions'][j] += int(len(my_seq2) / 2)
            else:
                p = df2.loc[i, 'Positions'][j]
                df2.loc[i, 'Positions'][j] = p - (int(len(my_seq2) / 2) + 1)
        j += 1
    df2.loc[i, 'Positions'].sort()
    dist = []
    distances(df2.loc[i, 'Positions'], dist)
    df2.loc[i, 'Distances'] = dist
    i += 1

# combine the 2 dataframes
frames = [df, df2]
result = pd.concat(frames)
df = result
df = df.reset_index(drop=True)


# change the lists to tuples, so identical rows can be deleted
i = 0
while i < len(df):
    pos = df.loc[i, 'Positions']
    dist = df.loc[i, 'Distances']
    tuple_pos = tuple(pos)
    tuple_dist = tuple(dist)
    df.loc[i, 'Positions'] = tuple_pos
    df.loc[i, 'Distances'] = tuple_dist
    i += 1

# deletes identical rows
df.drop_duplicates(inplace=True)
df = df.reset_index(drop=True)


# Delete redundant repeats
print("deleting redundant repeats")
i = 0
while i < len(df):
    deleted = False
    # creates new df with only repeats that contain the sequence of the row
    similar = df[df['Sequence'].str.contains(df.loc[i, 'Sequence'])]
    # if similar contains a repeat other than the repeat we're looking at
    if len(similar) > 1:
        # resets the index values of the similar df
        similar = similar.reset_index(drop=True)

        x = 0
        while x < len(similar) and not deleted:
            # check if this is the repeat we were looking at
            if similar.loc[x, 'Length'] > df.loc[i, 'Length']:
                # check if they have the same number of positions
                if similar.loc[x, 'Occurrences'] == df.loc[i, 'Occurrences']:
                    # delete the row we were looking at
                    df.drop([i], inplace=True)
                    df = df.reset_index(drop=True)
                    i -= 1
                    deleted = True
                    print("redundant")
            x += 1

    if not deleted:
        # check if 6 bp repeats overlap themselves
        if df.loc[i, 'Length'] == 6 and df.loc[i, "Occurrences"] == 2:
            if df.loc[i, 'Distances'][0] < 6:
                df.drop([i], inplace=True)
                df = df.reset_index(drop=True)
                i -= 1
                print("overlapping")
    print("rr: i = " + str(i) + " out of: " + str(len(df)-1))
    i += 1


# Sorts table by length and then sequence
df.sort_values(by=['Length', 'Sequence'], ascending=[False, True], inplace=True)
df = df.reset_index(drop=True)

# consolidate information for a sequence into 1 row
i = 0
print("consolidating sequences")
while i < len(df) - 1:
    if df.loc[i, 'Sequence'] == df.loc[i + 1, 'Sequence']:
        # updates positions
        pos1 = list(df.loc[i, 'Positions'])
        pos2 = list(df.loc[i + 1, 'Positions'])
        pos = pos1 + pos2
        pos = list(dict.fromkeys(pos))
        pos.sort()
        df.loc[i, 'Positions'] = tuple(pos)
        # updates distances
        dist = []
        distances(list(pos), dist)
        df.loc[i, 'Distances'] = tuple(dist)
        df.drop(i+1, inplace=True)
        # updates occurrences
        df.loc[i, 'Occurrences'] = len(pos)
        i += 1
    i += 1
df = df.reset_index(drop=True)


print(my_seq2)
print("\n")
print("final output table:")
print(df.to_string(index=False))
print("There are " + str(len(df)) + " repeated sequences")

# df.to_csv("pSL9 fasta scan.csv", index=False)
