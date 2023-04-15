import numpy as np
import pandas as pd
from SERGIO.sergio import sergio
import time
# import pickle
df = pd.read_csv('./bMat_cID7.tab', sep='\t', header=None, index_col=None)
bMat = df.values
sim = sergio(number_genes=1550, number_bins = 6, number_sc = 800, noise_params = 0.2, decays=0.8, sampling_state=1, noise_type='dpd', dynamics=True, bifurcation_matrix= bMat)
sim.build_graph(input_file_taregts ='5000.gene.6ct.reg.dyn.1.txt', 
                input_file_regs='5000.gene.6ct.m.reg.dyn.1.txt', 
                shared_coop_state=2)

st = time.time()
sim.simulate_dynamics()
et = time.time()
elapsed_time = et - st
print('Execution time:', elapsed_time, 'seconds')

exprU, exprS = sim.getExpressions_dynamics()
exprU_clean = np.concatenate(exprU, axis = 1)
exprS_clean = np.concatenate(exprS, axis = 1)

np.savetxt('5000.gene.6ct.dyn.v1.U.clean.csv', exprU_clean, delimiter = ',', fmt='%1.3f')
np.savetxt('5000.gene.6ct.dyn.v1.S.clean.csv', exprS_clean, delimiter = ',', fmt='%1.3f')




