#!/usr/bin/env python
# coding: utf-8

# In[8]:


import os
currect_directory = os.getcwd()
print(currect_directory)


# In[9]:


import pandas as pd
import numpy as np


# In[10]:


data = pd.read_csv('DecayTimecourse.txt', sep = '\t')


# In[11]:


print (data)


# In[12]:


timecourse_1 = pd.DataFrame(data[data.columns[1:10]])


# In[13]:


timecourse_1.head()


# In[14]:


timecourse_1.set_index(data['Time course #'], inplace = True)


# In[15]:


timecourse_1.columns = timecourse_1.iloc[0]


# In[16]:


timecourse_1.drop('YORF', axis = 0, inplace = True)


# In[17]:


timecourse_1.dropna(how = 'all', inplace = True)


# In[18]:


len(timecourse_1)


# In[19]:


log_timecourse_1 = np.log(timecourse_1)


# In[20]:


import statsmodels.api as sm


# In[21]:


results1 = {}
for i in range(len(log_timecourse_1.index)):
    if log_timecourse_1.iloc[i].isnull().any():
        x = log_timecourse_1.columns[~np.isnan(log_timecourse_1.iloc[i])]
        y = log_timecourse_1.iloc[i].dropna().values.reshape(-1, 1)
        x = sm.add_constant(x)
        model = sm.OLS(y, x).fit()
        slope = model.params[1]
        intercept = model.params[0]
        results1[log_timecourse_1.index[i]] = slope


# In[22]:


results1


# In[23]:


import numpy as np
# calculate the half-life for each row with missing values
half_life1 = {}
for i, slope in results1.items():
    half_life = np.log(2) / slope
    half_life1[i] = half_life

# print the half-life values
print(half_life1)


# In[24]:


df1 = pd.DataFrame.from_dict(half_life1, orient='index', columns=['half_life'])
# print the DataFrame
print(df1)


# In[25]:


timecourse_2 = pd.DataFrame(data[data.columns[10:19]])


# In[26]:


timecourse_2.head()


# In[27]:


timecourse_2.set_index(data['Time course #'], inplace = True)


# In[28]:


timecourse_2.columns = timecourse_2.iloc[0]


# In[29]:


timecourse_2.drop('YORF', axis = 0, inplace = True)


# In[30]:


timecourse_2.dropna(how = 'all', inplace = True)


# In[31]:


len(timecourse_2)


# In[32]:


log_timecourse_2 = np.log(timecourse_2)


# In[33]:


import statsmodels.api as sm


# In[34]:


results2 = {}
for i in range(len(log_timecourse_2.index)):
    if log_timecourse_2.iloc[i].isnull().any():
        x = log_timecourse_2.columns[~np.isnan(log_timecourse_2.iloc[i])]
        y = log_timecourse_2.iloc[i].dropna().values.reshape(-1, 1)
        x = sm.add_constant(x)
        model = sm.OLS(y, x).fit()
        slope = model.params[1]
        intercept = model.params[0]
        results2[log_timecourse_2.index[i]] = slope


# In[35]:


results2


# In[36]:


import numpy as np
# calculate the half-life for each row with missing values
half_life2 = {}
for i, slope in results2.items():
    half_life = np.log(2) / slope
    half_life2[i] = half_life

# print the half-life values
print(half_life2)


# In[37]:


df2 = pd.DataFrame.from_dict(half_life2, orient='index', columns=['half_life'])
# print the DataFrame
print(df2)


# In[38]:


timecourse_3 = pd.DataFrame(data[data.columns[19:28]])


# In[39]:


timecourse_3.head()


# In[40]:


timecourse_3.set_index(data['Time course #'], inplace = True)


# In[41]:


timecourse_3.columns = timecourse_3.iloc[0]


# In[42]:


timecourse_3.drop('YORF', axis = 0, inplace = True)


# In[43]:


timecourse_3.dropna(how = 'all', inplace = True)


# In[44]:


len(timecourse_3)


# In[45]:


log_timecourse_3 = np.log(timecourse_3)


# In[46]:


import statsmodels.api as sm


# In[47]:


results3 = {}
for i in range(len(log_timecourse_3.index)):
    if log_timecourse_3.iloc[i].isnull().any():
        x = log_timecourse_3.columns[~np.isnan(log_timecourse_3.iloc[i])]
        y = log_timecourse_3.iloc[i].dropna().values.reshape(-1, 1)
        x = sm.add_constant(x)
        model = sm.OLS(y, x).fit()
        slope = model.params[1]
        intercept = model.params[0]
        results3[log_timecourse_3.index[i]] = slope


# In[48]:


results3


# In[49]:


import numpy as np
# calculate the half-life for each row with missing values
half_life3 = {}
for i, slope in results3.items():
    half_life = np.log(2) / slope
    half_life3[i] = half_life

# print the half-life values
print(half_life3)


# In[50]:


df3 = pd.DataFrame.from_dict(half_life3, orient='index', columns=['half_life'])
# print the DataFrame
print(df3)


# In[51]:


import numpy as np

# calculate the average half-life for the three replicates
average_half_life = np.mean([df1, df2, df3])

# print the average half-life value
print('The average half-life is:', average_half_life)


# In[52]:


df4 = pd.DataFrame(average_half_life)
print(df4)


# In[60]:


# sort the genes based on their average half-life in ascending order
df_sorted = df4.sort_values('half_life')
# drop NaN values from the sorted dataframe
df_sorted = df4


# In[61]:


# print the sorted dataframe
print(df_sorted)


# In[62]:


# Calculate the number of rows that correspond to the top and bottom 10%
num_rows = len(df_sorted)
top_rows = int(num_rows * 0.1)
bottom_rows = int(num_rows * 0.1)

# Get the top and bottom 10% of genes
top_genes = df_sorted.head(top_rows)
bottom_genes = df_sorted.tail(bottom_rows)
print(top_genes)
print(bottom_genes)


# In[63]:


# save the top and bottom genes to text files
top_genes.to_csv('top_genes.csv')
bottom_genes.to_csv('bottom_genes.csv')

# download the files
from IPython.display import FileLink
display(FileLink('top_genes.csv'))
display(FileLink('bottom_genes.csv'))


# In[ ]:




