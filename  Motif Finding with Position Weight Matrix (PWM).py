#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
currect_directory = os.getcwd()
print(currect_directory)


# In[2]:


import numpy as np
import pandas as pd
# Load the counts matrix from the text file
counts_matrix = pd.read_csv('argR-counts-matrix.txt',sep='\t',header=None)
print(counts_matrix)


# In[3]:


counts_matrix.set_index(counts_matrix.columns[0], inplace=True)


# In[4]:


counts_matrix = counts_matrix.drop(counts_matrix.columns[0], axis=1)


# In[5]:


sum_series = counts_matrix.iloc[:, 0:].sum(axis=0)


# In[6]:


# Add a pseudocount of +1 to each count in the counts matrix
pseudocount_matrix = counts_matrix + 1
# Compute the total number of sequences (including pseudocounts)
num_seqs = np.sum(pseudocount_matrix)
# Compute the frequency matrix F'(b, j)
f_prime_matrix = pseudocount_matrix / num_seqs
# Compute the background frequency (assumed to be 0.25 for all bases)
bg_freq = 0.25
# Compute the frequency matrix F(b, j) using the formula W(b, j) = log(F'(b, j) / bg_freq) to represent weight 
weight_matrix = np.log2(f_prime_matrix / bg_freq)


# In[7]:


# Print the results
print("Frequency matrix F(b, j):")
print(weight_matrix)


# In[8]:


import numpy as np
import pandas as pd
df = pd.read_csv('E_coli_K12_MG1655.400_50',sep='\\',header=None)
print(df)


# In[9]:


df.drop(columns=[2], inplace=True)
print(df)


# In[10]:


df.columns = ['Gene_ID', 'Sequence']


# In[11]:


motif_length = len(weight_matrix.columns)
# write loop to iterate over each row extract Gene_ID and Sequence and identify motif score based on the weight matrix
results = {}
for index, row in df.iterrows():
    Gene_ID = row['Gene_ID']
    Sequence = row['Sequence'].replace(' ','')
    final_score = -np.inf
    for j in range(len(Sequence) - motif_length + 1):
        subsequence = Sequence[j:j+motif_length]
        score = sum(weight_matrix.loc[base, i+2] for i, base in enumerate(subsequence))
        if score > final_score:
            final_score = score
    results[Gene_ID] = final_score


# In[12]:


# Sort the genes by their scores in descending order
sorted_results = sorted(results.items(), key=lambda x: x[1], reverse=True)

# Print the top 30 genes
print("Top 30 genes:")
for gene, score in sorted_results[:30]:
    print(f"Gene ID: {gene}, Score: {score}")


# In[ ]:




