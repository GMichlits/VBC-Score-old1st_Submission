__author__ = 'georg.michlits'

import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split



#model_name = '3IMPsc_medm1_rel_v2'
#### add all improved tracr version screens
'''
Screen_data_namelist = ['/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/TKOv3_HAP1d18_CEG3_rel_LFC_d.sav']
Species = 'Hm'
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/Yusa/Exp_data_thirdCEG/Exp0third_rel_LFC_d.sav')
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/Yusa/Exp_data_thirdCEG/Exp1third_rel_LFC_d.sav')
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/Yusa/Exp_data_thirdCEG/Exp2third_rel_LFC_d.sav')
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/Yusa/Exp_data_thirdCEG/Exp3third_rel_LFC_d.sav')
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/Yusa/Exp_data_thirdCEG/Exp4third_rel_LFC_d.sav')
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/Yusa/Exp_data_thirdCEG/Exp5third_rel_LFC_d.sav')
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/Yusa/Exp_data_thirdCEG/Exp6third_rel_LFC_d.sav')
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/Yusa/Exp_data_thirdCEG/Exp7third_rel_LFC_d.sav')
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/Yusa/Exp_data_thirdCEG/Exp8third_rel_LFC_d.sav')
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/Yusa/Exp_data_thirdCEG/Exp9third_rel_LFC_d.sav')
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/Yusa/Exp_data_thirdCEG/Exp10third_rel_LFC_d.sav')
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/Yusa/Exp_data_thirdCEG/Exp11third_rel_LFC_d.sav')
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/Yusa/Exp_data_thirdCEG/Exp12third_rel_LFC_d.sav')
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Brunello_CEG3_rel_LFC_d.sav')
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/TKOv1_HTC116_CEG3_rel_LFC_d.sav')
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Wang15_Raji_Exp2CEG3_rel_LFC_d.sav')
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Wang15_K562_Exp3CEG3_rel_LFC_d.sav')
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Wang17_MOLM13_Exp136CEG3_rel_LFC_d.sav')
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Wang17_PL21_Exp128CEG3_rel_LFC_d.sav')
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/GeckoV2_Exp_12_HT29_LARGE_INTESTINE_CEG3_rel_LFC_d.sav')
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/GeckoV2_Exp_13_K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE_CEG3_rel_LFC_d.sav')
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Gecko_v2_Exp_25_PC3_PROSTATE_CEG3_rel_LFC_d.sav')
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Avana_Karpas_genCR_Exp285_CEG3_rel_LFC_d.sav')
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Avana_HCC1143_genCR_Exp227CEG3_rel_LFC_d.sav')
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Avana_NCIH1299_genCR_Exp347CEG3_rel_LFC_d.sav')
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Avana_RKO_genCR_Exp145CEG3_rel_LFC_d.sav')



Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/MIAPACA2_CEG3_rel_LFC_d.sav')
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/RKO_CEG3_rel_LFC_d.sav')
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/KBM7_CEG3_rel_LFC_d.sav')
'''
Screen_data_namelist = []
Species = 'Ms'
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/CrUMI_2n_CEG3_rel_LFC_d.sav')
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/CrUMI_n_CEG3_rel_LFC_d.sav')
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Zub_2n_d18_CEG3_rel_LFC_d.sav')
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Zub_n_d18_CEG3_rel_LFC_d.sav')
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Zub_2n_d6_CEG3_rel_LFC_d.sav')
Screen_data_namelist.append('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Zub_n_d6_CEG3_rel_LFC_d.sav')


#print(model_name + ' ' + data_name)

##### HartMoffat TKOv1
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/TKOv1_HTC116_CEG3_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/TKOv1_HTC116_abs_LFC_d.sav', 'rb'))
##### HartMoffat TKOv3
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/TKOv3_HAP1d18_CEG3_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/TKOv3_HAP1d18_abs_d.sav', 'rb'))
##### Wang15_Raji
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Wang15_Raji_Exp2CEG3_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Wang15_Raji_Exp2abs_LFC_d.sav', 'rb'))
##### Wang15_K562
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Wang15_K562_Exp3CEG3_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Wang15_K562_Exp3abs_LFC_d.sav', 'rb'))
##### Wang17_MOLM13
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Wang17_MOLM13_Exp136CEG3_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Wang17_MOLM13_Exp136abs_LFC_d.sav', 'rb'))
#####Wang17_PL21
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Wang17_PL21_Exp128CEG3_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Wang17_PL21_Exp128abs_LFC_d.sav', 'rb'))
##### Geckov2_HT29
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/GeckoV2_Exp_12_HT29_LARGE_INTESTINE_CEG3_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/GeckoV2_Exp_12_HT29_LARGE_INTESTINE_abs_d.sav', 'rb'))
##### Geckov2_K562
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/GeckoV2_Exp_13_K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE_CEG3_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/GeckoV2_Exp_13_K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE_abs_d.sav', 'rb'))
##### Geckov2 PC3
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Gecko_v2_Exp_25_PC3_PROSTATE_CEG3_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Gecko_v2_Exp_25_PC3_PROSTATE_abs_d.sav', 'rb'))

##### Yusa HL60 Exp3
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Yusa/Exp_data_thirdCEG/Exp3third_abs_LFC.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Yusa/Exp_data/Exp3abs_LFC_d.sav', 'rb'))
##### Yusa HT1080 Exp4
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Yusa/Exp_data_thirdCEG/Exp4third_abs_LFC.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Yusa/Exp_data/Exp4abs_LFC_d.sav', 'rb'))
##### Avana Karpas
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Avana_Karpas_genCR_Exp285_CEG3_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Avana_Karpas_genCR_Exp285_abs_LFC_d.sav', 'rb'))
##### Avana HCC1143 - Exp83
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Avana_HCC1143_genCR_Exp227CEG3_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Avana_HCC1143_genCR_Exp227abs_LFC_d.sav', 'rb'))
#### Brunello DATA LFC-0.5
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Brunello/Br_medm0.5_abs_LFC_d.sav', 'rb'))
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Brunello/Br_medm1_abs_LFC_d.sav', 'rb'))
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Brunello/Br_medm1.5_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Brunello/Brunello_abs_LFC_d.sav', 'rb'))
#### BRunello CEG3
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Brunello_CEG3_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Brunello_abs_LFC_d.sav', 'rb'))
#### IMP KBM7 DATA
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/KBM7_medm1_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/KBM7_abs_LFC_d.sav', 'rb'))
##### IMP RKO
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/RKO_medm1_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/RKO_abs_LFC_d.sav', 'rb'))
##### IMP MIApaca2
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/MIAPACA2_medm1_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/MIAPACA2_abs_LFC_d.sav', 'rb'))
##### CrUMI mESC Data_2n
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/CrUMI_2n_CEG3_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/CrUMI_2n_abs_LFC_d.sav', 'rb'))
##### CrUMI mESC Data_n
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/CrUMI_n_CEG3_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/CrUMI_n_abs_LFC_d.sav', 'rb'))
##### IMP mESC Data_2n
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Zub_2n_d18_CEG3_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Zub_2n_d18_abs_LFC_d.sav', 'rb'))
##### IMP mESC Data_n
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Zub_n_d18_CEG3_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Zub_n_d18_abs_LFC_d.sav', 'rb'))



#print('load model_name')
#model_T1_AAw5_6prop_abs_LFC = pickle.load(open('/Users/georg.michlits/Desktop/Models/VBC_score/model_Nov_rel.sav','rb'))
# gn30 Ms is gn63 Hm
#print('load Doench_with AApositional information')
#D16_AA_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/D16_AA_d.sav', 'rb'))
#print('load Doench_without AApositional information')
#D16_woAA_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/D16_woAA_d.sav', 'rb'))
#print('load inDelphi')
#inDelphi = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/inDelphi_infr_d.sav', 'rb'))

#inDelphi_fr_d = {}
#for pos in inDelphi:
#    inDelphi_fr_d[pos] = 1-inDelphi[pos]

#VBC_score_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/VBC_score_3IMPsc_CEG3_rel_v2_d.sav', 'rb'))
#VBC_score_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/VBC_score_KBM7_medm1_rel_v2_d.sav', 'rb'))
#VBC_score_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/VBC_score_KBM7_CEG3_rel_v2_d.sav', 'rb'))
#VBC_score_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/VBC_score_Nov_rel_v2_d.sav', 'rb'))
#VBC_score_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/VBC_score_3Nov_rel_v2_d.sav', 'rb'))
#VBC_score_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/VBC_score_3IMPsc_medm1_abs_v2_d.sav', 'rb'))
#VBC_score_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/VBC_score_3IMPsc_CEG3_abs_v2_d.sav', 'rb'))
#VBC_score_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/VBC_score_KBM7_medm1_abs_v2_d.sav', 'rb'))
#VBC_score_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/VBC_score_KBM7_CEG3_abs_v2_d.sav', 'rb'))



###################  Load basic information required
#print('load CEG')
#CEG_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/' + Species + '/property_sav/CEG_pos_d.sav', 'rb'))
print('load Gene_d')
Gene_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/' + Species + '/property_sav/Gene_d.sav', 'rb'))
Gene_direction_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/' + Species + '/property_sav/Gene_direction_d.sav', 'rb'))
print('load AA_residue_d')
#AA_residue_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/' + Species + '/property_sav/AA_residue_d+2.sav', 'rb'))
#AA_residue_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/' + Species + '/property_sav/AA_residue_d+1.sav', 'rb'))
#AA_residue_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/' + Species + '/property_sav/AA_residue_dv2.sav', 'rb'))
#AA_residue_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/' + Species + '/property_sav/AA_residue_d-1.sav', 'rb'))
AA_residue_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/' + Species + '/property_sav/AA_residue_d-2.sav', 'rb'))
###########################################################################

all_AA = 'RHKDESTNQGAVILMFYWCP'

def extract_bestworst_pos(Screen_data_rel_d, Gene_d):
    Struc_dict = {}
    best_worst = {}
    for pos in Screen_data_rel_d:
        if pos in Gene_d:
            rel_LFC = Screen_data_rel_d[pos]
            gene = Gene_d[pos]
            if not gene in Struc_dict:
                Struc_dict[gene] = {'pos': [], 'rel_LFC': []}
            Struc_dict[gene]['pos'].append(pos)
            Struc_dict[gene]['rel_LFC'].append(rel_LFC)
    for gene in Struc_dict:
        gene_pos_list = Struc_dict[gene]['pos']
        index_max = Struc_dict[gene]['rel_LFC'].index(np.max(Struc_dict[gene]['rel_LFC']))
        index_min = Struc_dict[gene]['rel_LFC'].index(np.min(Struc_dict[gene]['rel_LFC']))
        best_pos = gene_pos_list[index_max]
        worst_pos = gene_pos_list[index_min]
        best_worst[gene] = {'b':best_pos, 'w':worst_pos}
    return (best_worst)

score_matrix_output = open(Species + '_AAlist_outm2(S).txt','w')
#kcluster_matrix_output = open(Species + '_into_kcluster_out.txt','w')
#kcluster_matrix_output.write('Screen_Name\tdata')

score_matrix_output.write('screen')
timer = 0
for Screen_data_name in Screen_data_namelist:
    print(Screen_data_name + 'load')
    Screen_data_rel = pickle.load(open(Screen_data_name,'rb'))
    print(Screen_data_name + 'processed')

    AA_total = {}
    for letter in all_AA:
        AA_total[letter] = 0

    for pos in Screen_data_rel:
        guide_orient = pos.split('(')[1][0]
        if pos in Gene_direction_d:
            gene_orient = Gene_direction_d[pos]
            if pos in AA_residue_d and guide_orient == gene_orient:
                AA = AA_residue_d[pos]
                if AA not in AA_total:
                    AA_total[AA] = 0
                AA_total[AA] += 1

    AA_best_dict = {}
    AA_worst_dict = {}
    AA_score = {}

    for letter in all_AA:
        AA_best_dict[letter] = 0
        AA_worst_dict[letter] = 0
        AA_score[letter] = 0

    best_worst = extract_bestworst_pos(Screen_data_rel, Gene_d)
    for gene in best_worst:
        best_pos = best_worst[gene]['b']
        worst_pos = best_worst[gene]['w']
        guide_orient = best_pos.split('(')[1][0]
        if best_pos in Gene_direction_d:
            gene_orient = Gene_direction_d[best_pos]
            if best_pos in AA_residue_d and guide_orient == gene_orient:
                AA_best = AA_residue_d[best_pos]
                if AA_best not in AA_best_dict:
                    AA_best_dict[AA_best] = 0
                if AA_best not in AA_score:
                    AA_score[AA_best] = 0
                AA_best_dict[AA_best] += 1
                AA_score[AA_best] += 1
        guide_orient = worst_pos.split('(')[1][0]
        if worst_pos in Gene_direction_d:
            gene_orient = Gene_direction_d[worst_pos]
            if worst_pos in AA_residue_d and guide_orient == gene_orient:
                AA_worst = AA_residue_d[worst_pos]
                if AA_worst not in AA_worst_dict:
                    AA_worst_dict[AA_worst] = 0
                if AA_worst not in AA_score:
                    AA_score[AA_worst] = 0
                AA_worst_dict[AA_worst] += 1
                AA_score[AA_worst] += -1

    if timer == 0:
        for letter in all_AA:
            score_matrix_output.write('\t' + 'total_' + letter)
        for letter in all_AA:
            score_matrix_output.write('\t' + 'best_' + letter)
        for letter in all_AA:
            score_matrix_output.write('\t' + 'worst_' + letter)
        for letter in all_AA:
            score_matrix_output.write('\t' + 'score_' + letter)

    timer += 1
    score_matrix_output.write('\n' + Screen_data_name)
    for AA in all_AA:
        score_matrix_output.write('\t' + str(AA_total[AA]))
    for AA in all_AA:
        score_matrix_output.write('\t' + str(AA_best_dict[AA]))
    for AA in all_AA:
        score_matrix_output.write('\t' + str(AA_worst_dict[AA]))
    for AA in all_AA:
        score_matrix_output.write('\t' + str(AA_score[AA]))