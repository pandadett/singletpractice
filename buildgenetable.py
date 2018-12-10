import pandas as pd
import numpy as np
#Make sure to update the location - check for updates at biomart ensembl
gene_counts = pd.read_csv("/oak/stanford/groups/quake/nuttadap/learning_singlet/mart_export_new.txt",sep="\t")
#Drop duplicates of gene stable ID and only take the first
gene_counts_unique = gene_counts.drop_duplicates(subset='Gene stable ID',keep='first')
#Take only two components to create unique gene table
genetable = gene_counts_unique[['Gene stable ID','Gene name']]
gene_counts_unique_list = gene_counts_unique['Gene stable ID'].values
#Set index as 'Gene stable ID' 

#Import count table (this is from bag of stars)
dengue_counts = pd.read_csv("/oak/stanford/groups/quake/nuttadap/learning_singlet/elifedata/counts_dengue.tsv",sep='\t',index_col=0)
dcounts_gene = dengue_counts.index.values

# Develop a list of common genes:
common = []
for i in dcounts_gene:
        if i in gene_counts_unique_list:
               common.append(i)

# List of common genes + ERCC + NIST + special features
dcounts_common = dengue_counts[dengue_counts.index.isin(common)]+dengue_counts[dengue_counts.index.str.startswith("ER")]+dengue_counts[dengue_counts.index.str.startswith("__")]+dengue_counts[dengue_counts.index.str.startswith("N")]

# Make a table with common gene stable ID and gene names
gene_common = genetable[genetable['Gene stable ID'].isin(common)].sort_values(ascending=True,by = 'Gene stable ID')

# Make two dataframes:
df1 = pd.DataFrame()
df1['Gene stable ID'] = dcounts_common.index
df1 = df1.set_index('Gene stable ID')
df2 = gene_common
df2 = df2.set_index('Gene stable ID')
gene_name_final = pd.concat([df1,df2], axis=1, sort=True)
dcounts_common.insert(loc=0, column = "Gene name", value = gene_name_final['Gene name'])
