__author__ = 'georg.michlits'

import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split


print('load dataset')
RES16_LFC_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Novartis_tiled/DLD_abs_LFC_d.sav', 'rb'))

Species = 'Hm'
print('load Features')
D16_AA_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/D16_AA_d.sav', 'rb'))
AA_cons_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA7_d.sav', 'rb'))
distance_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/distance_d.sav', 'rb'))
mod3_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/mod3_d.sav', 'rb'))
VBC_score_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Apr1_VBC_score_all_rel_v2_d.sav', 'rb'))
Pfam_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Pfam_all_d.sav', 'rb'))
inDelphi_infr_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/inDelphi_infr_d.sav', 'rb'))
AA_score_dp1 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_residue_d+1.sav', 'rb'))
print('load other useful dictionaries')
Gene_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/' + Species + '/property_sav/Gene_d.sav', 'rb'))
nt30_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/' + Species + '/property_sav/nt30_d.sav', 'rb'))

Screen_data_all_d_list = []
Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/TKOv3_HAP1d18_abs_d.sav', 'rb')))
Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Wang17_MOLM13_Exp136abs_LFC_d.sav', 'rb')))
Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Avana_Karpas_genCR_Exp285_abs_LFC_d.sav', 'rb')))
Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Brunello_abs_LFC_d.sav', 'rb')))
Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/KBM7_abs_LFC_d.sav', 'rb')))
Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Yusa/Exp_data/Exp3abs_LFC_d.sav', 'rb')))

###################
#
#  NOTES
#   All dictionaries with _d contain "NAME_Property"_d = {pos : score}
#   All scores are between 0 and 1
#   pos is a guide ID string: chr12:85822107-85822137(+)
#   Dictionaries AA_cons_d has option of choosing window size "NAME_Property"_d = {pos : {1: score,
#                                                                                         3: score,
#                                                                                         5: score,.....}
####################

Data_list_dict = {}
Structure_dict = {}

print('generatefeatures for model')
i = 0
guide_not_found = 0
data = []
for pos in RES16_LFC_d:
    if pos in VBC_score_d:
        if pos in inDelphi_infr_d:
            data_list = []
            #####################################################################################
            ##                    x Features for MODEL
            AA_cons = VBC_score_d[pos]
            data_list.append(AA_cons)
            D16_AA = D16_AA_d[pos]
            data_list.append(D16_AA)
            inDelphi_infr = inDelphi_infr_d[pos]
            data_list.append(1-inDelphi_infr)
            Distance = distance_d[pos]
            data_list.append(Distance)
            mod3 = mod3_d[pos]
            data_list.append(mod3)
            Pfam = Pfam_all_d[pos]
            data_list.append(Pfam)

            gene = Gene_d[pos]
            Data_list_dict[pos] = data_list
            if gene not in Structure_dict:
                Structure_dict[gene] = {}
            Structure_dict[gene][pos] = RES16_LFC_d[pos]
        else:
            guide_not_found += 1
    else:
        guide_not_found += 1
    i += 1

print(str(i) + ' guides in dataset')
print(str(guide_not_found) + ' guides skipped. missing transcript ID - no AA_cons scores')

Feature_name_list = []
Feature_name_list.append('AA_cons')
Feature_name_list.append('gRNA act')
Feature_name_list.append('In-Frames')
Feature_name_list.append('distance')
Feature_name_list.append('mod3')
Feature_name_list.append('Pfam')
Target_name = 'Abs_LFC'

data_ML = {'data':np.array(data),
                'feature_names':np.array(Feature_name_list)}

print('features')
print(data_ML['feature_names'].shape)

for gene in Structure_dict:
    #os.mkdir(gene+'_Rb')
    i = 0
    x_list = []
    y_list = []
    c_list = []
    print(gene)
    for pos in sorted(Structure_dict[gene]):
        i += 1
        x_list.append(distance_d[pos])
        y_list.append(RES16_LFC_d[pos])
        c_list.append(round(Data_list_dict[pos][0], 4))
    sctr = plt.scatter(x_list, y_list, c=c_list, cmap='RdYlBu_r', linewidths=2)
    #plt.colorbar(sctr)
    plt.title(gene)
    plt.xlabel('guides along cDNA sequence')
    plt.ylabel('Log 2 fold change')
    plt.savefig('VBC_score'+'/' + gene + '_' + 'RB_r.pdf')
    plt.close()
    plt.close()

VBC_top10 = {}
for gene in Structure_dict:
    VBC_scores = {}
    for pos in Structure_dict[gene]:
        VBC_scores[pos] = VBC_score_d[pos]
    i = 0
    for item in sorted(VBC_scores.items(), key=lambda t: t[1], reverse=True):
        pos = item[0]
        VBC_score = item[1]
        if i < 10:
            VBC_top10[pos] = VBC_score
        i += 1
        if i == 10:
            break

Screen_data_all_d_list.append(VBC_top10)

for gene in Structure_dict:
    #os.mkdir(gene+'_Rb')
    i = 0
    x_list = []
    y_list = []
    c_list = []
    s_list = []
    print(gene)
    c_pick_list = []
    x_pick_list = []
    y_pick_list = []
    s_pick_list = []
    for pos in sorted(Structure_dict[gene]):
        i += 1
        x_list.append(distance_d[pos])
        y_list.append(RES16_LFC_d[pos])
        c_list.append('silver')
        s_list.append(20)
        screens_inanalysis = ['TKOv3','Sabatini','AVANA','Brunello','VBC','Sanger','VBCscore']
        screens_color = ['olive','lawngreen','green','b','c','dodgerblue','magenta']
        found = 0
        for n in range(7):
            if pos in Screen_data_all_d_list[n]:
                found += 15
                x_pick_list.append(distance_d[pos])
                y_pick_list.append(RES16_LFC_d[pos])
                c_pick_list.append(screens_color[n])
                s_pick_list.append(75-found)
                if gene == 'CTCF':
                    print(screens_inanalysis[n])
                    print(pos)
                    print(nt30_d[pos])
    sctr = plt.scatter(x_list, y_list, c=c_list,s=s_list)
    sctr = plt.scatter(x_pick_list, y_pick_list, c=c_pick_list, s=s_pick_list)
    plt.title(gene)
    plt.xlabel('guides along cDNA sequence')
    plt.ylabel('Log 2 fold change')
    plt.savefig('Screens_pick'+'/' + gene + '_' + '.pdf')
    plt.close()
    plt.close()
