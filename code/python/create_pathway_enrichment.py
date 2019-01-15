import pandas as pd
from bioservices import *
import re
import sys
import seaborn as sns
import matplotlib.pyplot as plt
import itertools
import numpy as np
from collections import Counter


def split_and_count(ser, col_name, delimiter=" "):

    # connect all strings and take care of NAs
    temp_s = delimiter.join([si for si in ser.values if isinstance(si, (str, unicode))])

    # split string which prepares for counting occurrences
    temp_l = temp_s.split()

    return pd.DataFrame(Counter(temp_l), index=[col_name]).transpose()


def get_ko_info(kegg_con, ko_df):

    ko_pw_rea = ko_df.copy()

    for koi in ko_pw_rea.index:

        ko_entry = kegg_con.parse(kegg_con.get(koi))

        try:
            ko_pw_rea.set_value(koi, 'pathways', " ".join(ko_entry['PATHWAY'].keys()))
        except (KeyError, TypeError):
            ko_pw_rea.set_value(koi, 'pathways', np.nan)


    return ko_pw_rea

# the input is KO_binary
org_ko_df = pd.read_csv('../../path/ko_binary.csv', index_col=0)

s = KEGG()

ko_infos = get_ko_info(s, org_ko_df)

col_keep = [coli for coli in ko_infos.columns if coli != 'pathways']

pw_count = pd.DataFrame()
for coli in col_keep:
    pw_temp = ko_infos[ko_infos[coli] == 1]['pathways']
    col_temp = split_and_count(pw_temp, coli)
    pw_count = pd.concat([pw_count, col_temp], axis=1)

pw_count.fillna(0, inplace=True)

for pi in pw_count.index:
    try:
        pw_count.set_value(pi, 'number_kos', len(s.parse(s.get(pi))['ORTHOLOGY'].keys()))
    except (KeyError, TypeError):
        print pi
        try:
            modules = s.parse(s.get(pi))['MODULE'].keys()
            print modules
            kos_pathw = set()
            for modi in modules:
                defi = s.parse(s.get(modi))['DEFINITION']
                kos_pathw.update(set(re.findall(r'[K]\d{5}', defi)))
            pw_count.set_value(pi, 'number_kos', len(kos_pathw))
        except (KeyError, TypeError):
            print pi
            pw_count.set_value(pi, 'number_kos', np.nan)
# OUTPUT
pw_count.to_csv('../../path/pw_count_ko_numbers.csv', index=True)
sys.exit(0)

# 3) map all KOs of a pathway on a map: http://www.genome.jp/kegg-bin/show_pathway?ko01230+K00030
# 3a) or can this be hacked? http://www.genome.jp/kegg/tool/map_pathway.html
