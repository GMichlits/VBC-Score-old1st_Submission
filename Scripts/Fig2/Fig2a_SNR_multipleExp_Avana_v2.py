__author__ = 'georg.michlits'

import pickle
import matplotlib.pyplot as plt
import os
import numpy as np

Species = 'Hm'
# enter name that a new generated folder will have where data will be stored
new_data_folder_name = 'Avana_GenomeCrispr_SNR'
#os.mkdir(new_data_folder_name)
# enter URL of folder all .sav files will be loaded and analysed all plots stored in /plots data stored in 'NAME_AUC.txt'
data_folder = '/Users/georg.michlits/Desktop/Gen_data/GenomeCrispr/Exp_data_Avana/'

CEG_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/' + Species + '/property_sav/CEG_pos_d.sav', 'rb'))
nonEG_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/' + Species + '/property_sav/nonEG_d.sav', 'rb'))

def det_mean_CEG_depl(Screen_data, CEG_d):
    i = 0
    j = 0
    k = 0
    CEG_depletions = []
    nonEG_depletions = []
    for pos in Screen_data:
        i += 1
        if pos in CEG_d:
            j += 1
            CEG_depletions.append(Screen_data[pos])
        if pos in nonEG_d:
            k += 1
            nonEG_depletions.append(Screen_data[pos])
    return(i,j,np.median(CEG_depletions),np.std(nonEG_depletions))

Avana_names_dict = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/GenomeCrispr/shortName_to_longName_Avana_d.sav', 'rb'))

for shorty in Avana_names_dict:
    print(shorty)
    print(Avana_names_dict[shorty])


AUC_outfile = open(new_data_folder_name + 'meanCEGdepl_v2.txt','w')
AUC_outfile.write('total_guides\tCEG_guides\tmean_CEG')
for File_name in os.listdir(data_folder):
    if File_name.endswith('LFC_d.sav'):
        print('read_file ' + File_name)
        data_name = File_name
        Screen_data = pickle.load(open(data_folder + '/' + data_name,'rb'))
        i, j, median_CEG_depl, std_nonEG_depl = det_mean_CEG_depl(Screen_data, CEG_d)
        AUC_outfile.write('\n' + Avana_names_dict[File_name.split('abs_')[0]] + '\t' + str(i) + '\t' + str(j) + '\t' + str(median_CEG_depl) + '\t' + str(std_nonEG_depl))