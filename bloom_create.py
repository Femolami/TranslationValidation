#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 13:54:38 2020

@author: bkearney
"""

"""

@author: bkearney
email: bkearne5@uncc.edu

Purpose: improve from TranslationfromCSV code with new datasets and more ubiquitous code

"""

import numpy as np
import pandas as pd
import re
import os
import statistics 
from pathlib import Path
import os


# import allel

# Change directory to path with all reference/input files
    
gene = pd.read_table('gene2ensembl')
human = pd.read_table('Homo_sapiens.gene_info')
mart = pd.read_csv('mart_export.txt',delimiter=',')

# Create Aliases column (delimited Synonyms)
# aliases = human['Synonyms'].tolist()
# aliases_sep = []
# for i in range(len(aliases)):
#     my_list = aliases[i].split("|")
#     my_list.append(human['Symbol'].iloc[i])
#     aliases_sep.append(my_list)
# human['Aliases'] = aliases_sep

# Create geneId column from splicing/filtering dbXrefs
dbXrefs = human['dbXrefs'].tolist()
geneIds = []
for i in range(len(dbXrefs)):
    if 'Ensembl' in dbXrefs[i]:
        gene_id = dbXrefs[i].split("Ensembl:")[1]
    else:
        gene_id = ""
    geneIds.append(gene_id)
human['geneId'] = geneIds

human_modified = human.drop(columns=['GeneID','LocusTag','Modification_date','Feature_type'])
mart_modified = mart.drop(columns=['Gene start (bp)','Gene end (bp)','Transcript start (bp)','Transcript end (bp)'])

from bloom_filter import BloomFilter

def matching_check(name, bloom):
    if len(name)>0:
        if name in bloom:
            return True
                
        else:
            return False
    
    
def add_to_bloom(df, bloom_name):
    for row in df.columns:
        for ele in df[row]:
            if (matching_check(str(ele),bloom_name) is False):
                bloom_name.add(str(ele))
                

gene_bloom = BloomFilter(max_elements=100000000)
add_to_bloom(gene, gene_bloom)

human_bloom = BloomFilter(max_elements=10000000)
add_to_bloom(human_modified,human_bloom)

mart_bloom = BloomFilter(max_elements=10000000)
add_to_bloom(mart,mart_bloom)

blooms = [human_bloom, gene_bloom, mart_bloom]
bloom_names = ['human','gene','mart']



def bloom_translate(file):
    df = pd.read_excel(file)

    cols = df.columns
    
    list_res = []
    list_label = []
    counter = 0
    for ele in blooms:
        final_match = []
        for col_name in cols:
            match = 0
            miss = 0
            for i in range(len(df[col_name])):
                name = str(df[col_name][i])
                if matching_check(name,ele):
                    match=match+1
                else:
                    miss = miss+1
            final = match/(match+miss)
            final_match.append(final)
#        res = {cols[i]: final_match[i] for i in range(len(cols))} 
        res = {cols[i]: final_match[i] for i in range(len(cols))} 
        res = {key:val for key, val in res.items() if val > .8}
#        list_res.append(list(res.values()))
        list_label.append(len(list(res.keys())))
        counter = counter+1
    # List what it matched, source file name and column, top 3 hits
#    return list_res,list_label
    return list_label

def process_files():
    rootdir = Path('/Users/bkearney/Documents/R_projects/pmc_xls')
#    file_list = [f for f in rootdir.glob('**/*') if f.is_file()]
#    print(type(file_list[0]))    

    total_list = []
    counter = 0
    for subdir, dirs, files in os.walk(rootdir):
        for file in files:
            try:
                file_name = os.path.join(subdir, file)
                translated_file = bloom_translate(file_name)
#                print(file_name,translated_file[1])
                total_list.append(translated_file)
                print(translated_file)
                print("Counter: ",counter)
                counter = counter+1
            except Exception:
                pass
            
            

            
#    pmc507897 = bloom_translate("/Users/bkearney/Documents/R_projects/pmc_xls/PMC507879/gb-2004-5-8-r54-s1.xls")
    df_final = pd.DataFrame(total_list, columns = ['gene','human','mart'])
    return df_final
#drug = bloom_translate('drug_supplement.xlsx')


df90 = process_files()
print(df90)
    # List what it matched, source file name and column, top 3 hits

# Manually confirm a few/subset files and compare contents to what got matched
# Estimate how accurate it is overall