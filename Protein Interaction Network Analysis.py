#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
currect_directory = os.getcwd()
print(currect_directory)


# In[2]:


import networkx as nx
import pandas as pd
# Load the edge list from the text file
network_df = pd.read_csv('Human-PPI.txt', sep='\s+')


# In[3]:


network_df.head()


# In[4]:


import networkx as nx
# create an empty undirected graph
graph = nx.Graph()


# In[5]:


# Iterate over the rows of the pandas dataframe
for _, row in network_df.iterrows():
    # add the two proteins as nodes to the graph
    graph.add_node(row['OFFICIAL_SYMBOL_A'])
    graph.add_node(row['OFFICIAL_SYMBOL_B'])
    # add the interaction as an edge between the two nodes
    graph.add_edge(row['OFFICIAL_SYMBOL_A'], row['OFFICIAL_SYMBOL_B'])


# In[6]:


print('Number of nodes in the graph:', graph.number_of_nodes())


# In[7]:


print('Number of edges in the graph:', graph.number_of_edges())


# In[8]:


# calculate the degree and clustering coefficient for each node
degrees = nx.degree(graph)
clustering = nx.clustering(graph)


# In[9]:


# calculate the average clustering coefficient of the network
avg_clustering = nx.average_clustering(graph)


# In[10]:


print("Average clustering coefficient:", avg_clustering)


# In[11]:


# get the degree sequence of the graph
degree_sequence = sorted([d for n, d in degrees], reverse=True)


# In[12]:


import numpy as np
import matplotlib.pyplot as plt

# Get the degree of each node
degrees = dict(graph.degree())

# Get the degree values and their frequency
degree_values = np.array(list(degrees.values()))
degree_freq = np.bincount(degree_values)

# Plot the degree distribution on a log-log scale
plt.loglog(np.arange(len(degree_freq)), degree_freq, 'o')
plt.xlabel('Degree')
plt.ylabel('Frequency')
plt.title('Degree Distribution')
plt.show()


# In[13]:


import networkx as nx
import pandas as pd
from scipy.stats import wilcoxon
from scipy.stats import ranksums


# In[14]:


# read the file into a pandas DataFrame
proteins1 = pd.read_csv('Protein-list1.txt')

# rename the 'LAS1L' column to 'protein A'
proteins1 = proteins1.rename(columns={'LAS1L': 'protein A'})

# insert a new row with the value 'LAS1L' in the first column
proteins1.loc[-1] = ['LAS1L'] + [''] * (len(proteins1.columns) - 1)

# reset the index to start from 0
proteins1 = proteins1.reset_index(drop=True)

# print the modified DataFrame to verify the change
print(proteins1)


# In[15]:


# read the file into a pandas DataFrame
proteins2 = pd.read_csv('Protein-list2.txt')

# rename the 'GNL3L' column to 'protein B'
proteins2 = proteins2.rename(columns={'GNL3L': 'protein B'})

# insert a new row with the value 'GNL3L' in the first column
proteins2.loc[-1] = ['GNL3L'] + [''] * (len(proteins2.columns) - 1)

# reset the index to start from 0
proteins2 = proteins2.reset_index(drop=True)

# print the modified DataFrame to verify the change
print(proteins2)


# In[16]:


# create subgraphs containing only proteins in protein list 1 and protein list 2
subgraph1 = nx.Graph()
subgraph2 = nx.Graph()


# In[17]:


subgraph1_path_lengths = []
for i in range(len(proteins1)-1):
    for j in range(i+1, len(proteins1)):
        p1 = proteins1['protein A'][i]
        p2 = proteins1['protein A'][j]
        if p1 in list(graph.nodes()) and p2 in list(graph.nodes()):
            subgraph1_path_lengths.append(len(nx.shortest_path(graph,p1,p2)))
        else:
            continue


# In[18]:


subgraph2_path_lengths = []
for i in range(len(proteins2)-1):
    for j in range(i+1, len(proteins2)):
        p1 = proteins2['protein B'][i]
        p2 = proteins2['protein B'][j]
        if p1 in list(graph.nodes()) and p2 in list(graph.nodes()):
            subgraph2_path_lengths.append(len(nx.shortest_path(graph,p1,p2)))
        else:
            continue


# In[19]:


wilcoxon_test = ranksums(subgraph1_path_lengths, subgraph2_path_lengths)


# In[20]:


wilcoxon_test


# In[ ]:




