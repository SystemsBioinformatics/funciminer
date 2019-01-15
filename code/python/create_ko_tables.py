import pandas as pd
import os
import glob
import re
import sys
import seaborn as sns
import matplotlib.pyplot as plt


# path to KEGG functional annotation files
all_dir = glob.glob("../../path/*_KO.txt")

df_all_kos = pd.DataFrame()

# sys.exit(0)

for kof in all_dir:

    try:
        inp = open(kof, 'r')

    except IOError:
        print 'File does not exist'

    # read the KO table and assign names to columns
    df_temp = pd.read_csv(inp, names=['pid', 'ko'], sep='\t', header=None)
    # get rid of NAs
    df_temp.dropna(inplace=True)
    # count occurrences of KOs
    ser_counts = df_temp['ko'].value_counts()
    # create new data frame with counts
    df_trans = pd.DataFrame(ser_counts).rename(columns={'ko': kof.rsplit('/', 2)[-1][:-7]}).transpose()
    # concatenate dataframes; if KOs are not present, NAs will be introduced
    df_all_kos = pd.concat([df_all_kos, df_trans], axis=0)

# concat will lead to NAs if an organism does not have certain KOs; replace them by 0
df_all_kos.fillna(0, inplace=True)
#FIRST OUTPUT
df_all_kos.transpose().to_csv('../../path/ko_counts.csv', index=True)

# binary representation
df_kos_binary = df_all_kos.copy()
df_kos_binary[df_kos_binary > 0] = 1
df_kos_binary = df_kos_binary.T
# SECOND OUTPUT
df_kos_binary.to_csv('../../path/ko_binary.csv', index=True)


sys.exit(0)

# replace IDs by actual names
name_df = pd.read_excel('../../data/gb_fasta_files/kefir/Bacteria_isolates.xlsx')

for coli in list(df_kos_binary):
    # print 'column: '+coli
    for entryi in name_df['Colony no.']:
        # print 'entry: '+str(entryi)
        if coli.startswith(str(entryi)+'_'):
            # print 'entry string: '+str(entryi)+'_'
            df_kos_binary.rename(columns={coli: name_df.Species[name_df['Colony no.'] == entryi].values.tolist()[0]},
                                 inplace=True)

df_kos_binary.to_csv('../../data/ko_tables/kefir/ko_binary_actual_names.csv', index=True)

sys.exit(0)

