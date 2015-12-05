
# coding: utf-8

# In[151]:

get_ipython().magic('matplotlib inline')
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
pd.set_option('display.mpl_style', 'default') 
pd.set_option('display.width', 5000) 
pd.set_option('display.max_columns', 60)



# In[141]:

flat_coalescence_stats =     pd.read_csv        ('./Dropbox/Thesis/G-PhoCS/code/G-PhoCS/logs/sample-data.coalStats.tsv',         sep='\t')


# In[148]:

flat_coalescence_stats['P_Z_ϴM'] = flat_coalescence_stats['dataLogLikelihood'] - flat_coalescence_stats['genLogLikelihood']


# In[160]:

np.exp(-62043.527880)
#flat_coalescence_stats[['P_Z_ϴM']].applymap(np.exp)


# In[ ]:



