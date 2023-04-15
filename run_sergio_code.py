import numpy as np
import pandas as pd
from SERGIO.sergio import sergio
import time

sim = sergio(number_genes=5000, number_bins = 3, number_sc = 1500, noise_params = 1, decays=0.8, sampling_state=15, noise_type='dpd')
sim.build_graph(input_file_taregts ='5000.gene.3ct.reg.1.txt', 
                input_file_regs='5000.gene.3ct.m.reg.1.txt', 
                shared_coop_state=2)

st = time.time()
sim.simulate()
et = time.time()
elapsed_time = et - st
print('Execution time:', elapsed_time, 'seconds')


expr = sim.getExpressions()
expr_clean = np.concatenate(expr, axis = 1)
np.savetxt('5000.gene.3ct.v1.clean.csv', expr_clean, delimiter = ',', fmt='%1.3f')



