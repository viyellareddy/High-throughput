#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
currect_directory = os.getcwd()
print(currect_directory)


# In[2]:


import pandas as pd
# Load PTT data
data1 = pd.read_csv('E_coli_K12_MG1655.ptt', sep='\t', skiprows=[0,1])
print(data1)


# In[3]:


def predict_operons(data):
    operons = []
    current_operon = []
    for _, row in data.iterrows():
        if current_operon == []:
            current_operon.append(row)
        elif row['Strand'] == current_operon[-1]['Strand'] and int(row['Location'].split('..')[0]) - int(current_operon[-1]['Location'].split('..')[1]) < 50:
            current_operon.append(row)
        else:
            operons.append(current_operon)
            current_operon = [row]
    operons.append(current_operon)
    return operons


# In[4]:


# Predict operons
operons = predict_operons(data1)
for i, operon in enumerate(operons):
    print(f"Operon {i+1}: {', '.join(gene['Gene'] for gene in operon)}")


# In[5]:


data2 = pd.read_csv('B_subtilis_168.ptt', sep='\t', skiprows=[0,1])
print(data2)


# In[6]:


def predict_operons(data):
    operons = []
    current_operon = []
    for _, row in data.iterrows():
        if current_operon == []:
            current_operon.append(row)
        elif row['Strand'] == current_operon[-1]['Strand'] and int(row['Location'].split('..')[0]) - int(current_operon[-1]['Location'].split('..')[1]) < 50:
            current_operon.append(row)
        else:
            operons.append(current_operon)
            current_operon = [row]
    operons.append(current_operon)
    return operons


# In[7]:


# Predict operons
operons = predict_operons(data2)
for i, operon in enumerate(operons):
    print(f"Operon {i+1}: {', '.join(gene['Gene'] for gene in operon)}")


# In[8]:


data3 = pd.read_csv('Halobacterium_NRC1.ptt', sep='\t', skiprows=[0,1])
print(data3)


# In[9]:


def predict_operons(data):
    operons = []
    current_operon = []
    for _, row in data.iterrows():
        if current_operon == []:
            current_operon.append(row)
        elif row['Strand'] == current_operon[-1]['Strand'] and int(row['Location'].split('..')[0]) - int(current_operon[-1]['Location'].split('..')[1]) < 50:
            current_operon.append(row)
        else:
            operons.append(current_operon)
            current_operon = [row]
    operons.append(current_operon)
    return operons


# In[10]:


# Predict operons
operons = predict_operons(data3)
for i, operon in enumerate(operons):
    print(f"Operon {i+1}: {', '.join(gene['Gene'] for gene in operon)}")


# In[11]:


data4 = pd.read_csv('Synechocystis_PCC6803_uid159873.ptt', sep='\t', skiprows=[0,1])
print(data4)


# In[12]:


def predict_operons(data):
    operons = []
    current_operon = []
    for _, row in data.iterrows():
        if current_operon == []:
            current_operon.append(row)
        elif row['Strand'] == current_operon[-1]['Strand'] and int(row['Location'].split('..')[0]) - int(current_operon[-1]['Location'].split('..')[1]) < 50:
            current_operon.append(row)
        else:
            operons.append(current_operon)
            current_operon = [row]
    operons.append(current_operon)
    return operons


# In[13]:


# Predict operons
operons = predict_operons(data4)
for i, operon in enumerate(operons):
    print(f"Operon {i+1}: {', '.join(gene['Gene'] for gene in operon)}")


# In[14]:


import pandas as pd
# Load GFF data for Hoatzin crop microbiome
data = pd.read_csv('2088090036.gff', sep='\t', header=None, comment ='#')
data.columns = ['Contig', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes']
# Sort data by contig and start position
data = data.sort_values(['Contig', 'Start'])


# In[15]:


print(data)


# In[16]:


# Extract gene names from attributes column
data['Gene'] = data['Attributes'].apply(lambda x: x.split(';')[0].split('=')[1])

# Define function to predict operons
def predict_operons(data):
    operons = []
    current_operon = []
    for _, row in data.iterrows():
        if current_operon == []:
            current_operon.append(row)
        elif row['Strand'] == current_operon[-1]['Strand'] and int(row['Start']) - int(current_operon[-1]['End']) < 50:
            current_operon.append(row)
        else:
            operons.append(current_operon)
            current_operon = [row]
    operons.append(current_operon)
    return operons


# In[17]:


# Predict operons
operons = predict_operons(data)
# Print predicted operons
for i, operon in enumerate(operons):
    print(f"Operon {i+1}: {', '.join(row['Gene'] for row in operon)}")


# In[ ]:




