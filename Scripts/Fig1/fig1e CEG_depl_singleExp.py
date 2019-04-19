__author__ = 'georg.michlits'

import pickle
import matplotlib.pyplot as plt
import numpy as np

Species = 'Hm'
data_folder = ''    #enter URL of folder all .sav files will be loaded and analysed all plots stored in /plots data_name stored in 'NAME_AUC.txt'
data = 'Brunello'
' + Species + '
Screen_data = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Brunello_abs_LFC_d.sav','rb'))
CEG_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/' + Species + '/property_sav/CEG_pos_d.sav', 'rb'))

i = 0
j = 0
CEG_depletions = []
for pos in Screen_data:
    i += 1
    if pos in CEG_d:
        j += 1
        CEG_depletions.append(Screen_data[pos])

print(i)
print(j)
print(np.mean(CEG_depletions))

AUC_outfile = open(data + 'meanCEGdepl.txt','w')
AUC_outfile.write('total_guides\tCEG_guides\tmean_CEG')
AUC_outfile.write('\n' + str(i) + '\t' + str(j) + '\t' + str(np.mean(CEG_depletions)))
