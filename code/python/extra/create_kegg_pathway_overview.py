#!/usr/bin/env python

#create overview table over pathways in kegg (number of reactions, numberof kos, number od compounds, number of modules)

import os, sys
import xml.etree.ElementTree as ET
import requests
from bioservices import KEGG


def get_num_reactions(ko):
    reactions = s.link('reaction', 'map'+ko[2:]).split('\n')
    return len([str(r) for r in reactions if r != ''])
    
def get_num_kos(ko):
    kos = s.link('ko', 'map'+ko[2:]).split('\n')
    return len([str(k) for k in kos if k != ''])
    
def get_num_compounds(ko):
    compounds = s.link('compound', 'map'+ko[2:]).split('\n')
    return len([str(c) for c in compounds if c != ''])
    
def get_num_modules(ko):
    modules = s.link('module', 'map'+ko[2:]).split('\n')
    return len([str(m) for m in modules if m != ''])
    

def get_all_pathway_ids_from_kegg(id_out_file):
    """ gets the pathway ids from kegg and counts the number of reactions, compounds, kos and modules for each pathway """
    with open(id_out_file, 'w') as out:
        kegg_pathway_list = requests.get('http://rest.kegg.jp/list/pathway')
        out.write('id\tname\t#compounds\t#reactions\t#kos\t#modules\n')     
        try:   
            for line in kegg_pathway_list.content.split('\n'):
                pwid, pwname = line.split('\t')
                pwid = 'ko' + pwid[8:]
                rcount = get_num_reactions(pwid)
                ccount = get_num_compounds(pwid)
                kcount = get_num_kos(pwid)
                mcount = get_num_modules(pwid)
                txt = '%s\t%s\t%d\t%d\t%d\t%d\n' % (pwid, pwname, ccount, rcount, kcount, mcount)
                out.write(txt)
        except ValueError:
            pass
            

s = KEGG()
id_out_file = 'kegg_pathways_overview.tab'
if len(sys.argv) > 1:
    id_out_file = sys.argv[1]
get_all_pathway_ids_from_kegg(id_out_file)

