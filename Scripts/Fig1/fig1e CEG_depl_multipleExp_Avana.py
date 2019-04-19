__author__ = 'georg.michlits'

import pickle
import matplotlib.pyplot as plt
import os
import numpy as np

Species = 'Hm'
# enter name that a new generated folder will have where data will be stored
new_data_folder_name = 'Avana_Hm'
os.mkdir(new_data_folder_name)
# enter URL of folder all .sav files will be loaded and analysed all plots stored in /plots data stored in 'NAME_AUC.txt'
data_folder = '/Users/georg.michlits/Desktop/Gen_data/Avana_cell_lines/all_LFCv1/'

CEG_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/' + Species + '/property_sav/CEG_pos_d.sav', 'rb'))
nonEG_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/' + Species + '/property_sav/nonEG_d.sav', 'rb'))

def det_mean_CEG_depl(Screen_data, CEG_d):
    i = 0
    j = 0
    CEG_depletions = []
    for pos in Screen_data:
        i += 1
        if pos in CEG_d:
            j += 1
            CEG_depletions.append(Screen_data[pos])
    return(i,j,np.mean(CEG_depletions))

AUC_outfile = open(new_data_folder_name + 'meanCEGdepl.txt','w')
AUC_outfile.write('total_guides\tCEG_guides\tmean_CEG')
for File_name in os.listdir(data_folder):
    if File_name.endswith('LFC_d.sav'):
        print('read_file ' + File_name)
        data_name = File_name
        Screen_data = pickle.load(open(data_folder + '/' + data_name,'rb'))
        i, j, mean_CEG_depl = det_mean_CEG_depl(Screen_data, CEG_d)
        AUC_outfile.write('\n' + str(i) + '\t' + str(j) + '\t' + str(mean_CEG_depl))