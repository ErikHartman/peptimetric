import pandas as pd

def read_file(file): # reads a file and outputs a dataframe
    df = pd.read_excel(file)
    return df

def drop_zeros(df, colname): #drops 0s in dataframe for given column
    df = df[df[colname] != 0]
    return df

def amino_acid_frequency(list): #gets the frequency for amino acids
    letters = {
        'A': 0,
        'G': 0,
        'V': 0,
        'L': 0,
        'I': 0,
        'P': 0,
        'F': 0,
        'W': 0,
        'M': 0,
        'S': 0,
        'T': 0,
        'C': 0,
        'Y': 0,
        'N': 0,
        'Q': 0,
        'K': 0,
        'R': 0,
        'H': 0,
        'D': 0,
        'E': 0
    }
    for sequence in list:
        for letter in sequence:
            letters[letter] += 1
    return letters

def group(list): #Groups amino acids for sequences in a list. Returns grouped
    grouped=[]
    nonpolar=['G','A','V','L','I','P','F','W','M']
    polar=['S','T','C','Y','N','Q']
    basic=['K','R','H']
    acidic=['D','E']
    for sequence in list:
        new_item=''
        for letter in sequence:
            if letter in nonpolar:
                new_item+='N'
            if letter in polar:
                new_item+='P'
            if letter in basic:
                new_item+='B'
            if letter in acidic:
                new_item+='A'
        grouped.append(new_item)
    return grouped


def print_my_list(list):
    print(list)