__author__ = 'georg.michlits'

import numpy as np
import matplotlib.pyplot as plt
import pickle

print('load data')
Screen = 'MIAPACA2'
Species = 'Hm'
LFC_data = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/'+Screen+'_abs_LFC_d.sav','rb'))
ppos = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/'+Screen+'_ppos_d.sav', 'rb'))
pneg = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/'+Screen+'_pneg_d.sav', 'rb'))
CEG_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/CEG_pos_d.sav', 'rb'))
plasmid_count = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/Zub_hm_plasmid_d.sav', 'rb'))

x_ls = []
y_ls = []
x_CEG = []
y_CEG = []
for pos in LFC_data:
    LFC = LFC_data[pos]
    x_ls.append(LFC)
    y_ls.append(plasmid_count[pos])
    if pos in CEG_d:
        x_CEG.append(LFC)
        y_CEG.append(plasmid_count[pos])

print(np.array(x_ls))
print(np.array(y_ls))
plt.figure(figsize=(12,10))
plt.yscale('log')
plt.scatter(x_ls,y_ls, c='grey', s=30, alpha=0.3)
plt.scatter(x_CEG,y_CEG, c='red', s=30, alpha=0.3)
plt.xlim(-12,4)
plt.tick_params(labelsize=20)
#plt.show()
plt.savefig(Screen + '_vs_libReads.png')
