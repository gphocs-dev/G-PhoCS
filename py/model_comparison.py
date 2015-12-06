get_ipython().magic('matplotlib inline')
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
pd.set_option('display.mpl_style', 'default') 
pd.set_option('display.width', 5000) 
pd.set_option('display.max_columns', 60)



flat_coalescence_stats =     pd.read_csv        ('../logs/sample-data.flatStats.tsv',         sep='\t')




flat_coalescence_stats = flat_coalescence_stats[['logPrior','coalStatFlat', 'numCoalFlat', 'root_theta', 'genealogyLogLikelihood' ]]
flat_coalescence_stats.columns = ['logPrior','time_stats', 'num_coal', 'root_ϴ', 'P_Z_ϴM' ]

printFactor = 1000.0
flat_coalescence_stats['root_ϴ'] = flat_coalescence_stats[['root_ϴ']].apply(lambda x:x/printFactor)

flat_coalescence_stats[:5]




def P_Z_ϴM0(theta, num_coal, time_stats):
    result = num_coal*np.log(2.0/theta) -(time_stats/theta)
    return result
P_Z_ϴM0(0.00282063, 7000, 3.442481)




flat_coalescence_stats['root_ϴ'] = flat_coalescence_stats[['root_ϴ']].apply(lambda x:x/1000.0)
root_ϴ = flat_coalescence_stats['root_ϴ'] 
num_coal = flat_coalescence_stats['num_coal'] 
time_stats = flat_coalescence_stats['time_stats']

flat_coalescence_stats['P_Z_ϴM0'] = P_Z_ϴM0(root_ϴ, num_coal, time_stats)

flat_coalescence_stats.sort_values(by='P_Z_ϴM0', ascending=False)[:5]


# In[ ]:



