__author__ = 'georg.michlits'

import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
'3Nov_rel_v2'
'3IMPsc_medm1_rel_v2'
model_name_list = ['all_rel_v2']
#model_name = '3IMPsc_medm1_rel_v2'
Analysis_name = 'Apr1_Ms_scorercomp'

#print(model_name + ' ' + data_name)

##################################################################################################################################
#################################################### HUMAN
'''
Species = 'Hm'
Screen_name_list = []
Screen_data_d_list = []
Screen_data_all_d_list = []
Screen_name_list.append('HartMoffat TKOv1')
Screen_data_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/TKOv1_HTC116_CEG3_abs_LFC_d.sav', 'rb')))
Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/TKOv1_HTC116_abs_LFC_d.sav', 'rb')))
Screen_name_list.append('HartMoffat TKOv3')
Screen_data_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/TKOv3_HAP1d18_CEG3_abs_LFC_d.sav', 'rb')))
Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/TKOv3_HAP1d18_abs_d.sav', 'rb')))
Screen_name_list.append('Wang15_Raji')
Screen_data_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Wang15_Raji_Exp2CEG3_abs_LFC_d.sav', 'rb')))
Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Wang15_Raji_Exp2abs_LFC_d.sav', 'rb')))
Screen_name_list.append('Wang15_K562')
Screen_data_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Wang15_K562_Exp3CEG3_abs_LFC_d.sav', 'rb')))
Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Wang15_K562_Exp3abs_LFC_d.sav', 'rb')))
Screen_name_list.append('Wang17_MOLM13')
Screen_data_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Wang17_MOLM13_Exp136CEG3_abs_LFC_d.sav', 'rb')))
Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Wang17_MOLM13_Exp136abs_LFC_d.sav', 'rb')))
Screen_name_list.append('Wang17_PL21')
Screen_data_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Wang17_PL21_Exp128CEG3_abs_LFC_d.sav', 'rb')))
Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Wang17_PL21_Exp128abs_LFC_d.sav', 'rb')))
Screen_name_list.append('Geckov2_HT29')
Screen_data_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/GeckoV2_Exp_12_HT29_LARGE_INTESTINE_CEG3_abs_LFC_d.sav', 'rb')))
Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/GeckoV2_Exp_12_HT29_LARGE_INTESTINE_abs_d.sav', 'rb')))
Screen_name_list.append('Geckov2_K562')
Screen_data_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/GeckoV2_Exp_13_K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE_CEG3_abs_LFC_d.sav', 'rb')))
Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/GeckoV2_Exp_13_K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE_abs_d.sav', 'rb')))
Screen_name_list.append('Geckov2 PC3')
Screen_data_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Gecko_v2_Exp_25_PC3_PROSTATE_CEG3_abs_LFC_d.sav', 'rb')))
Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Gecko_v2_Exp_25_PC3_PROSTATE_abs_d.sav', 'rb')))
Screen_name_list.append('Yusa HL60 Exp3')
Screen_data_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Yusa/Exp_data_thirdCEG/Exp3third_abs_LFC.sav', 'rb')))
Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Yusa/Exp_data/Exp3abs_LFC_d.sav', 'rb')))
Screen_name_list.append('Yusa HT1080 Exp4')
Screen_data_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Yusa/Exp_data_thirdCEG/Exp4third_abs_LFC.sav', 'rb')))
Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Yusa/Exp_data/Exp4abs_LFC_d.sav', 'rb')))
Screen_name_list.append('Avana Karpas')
Screen_data_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Avana_Karpas_genCR_Exp285_CEG3_abs_LFC_d.sav', 'rb')))
Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Avana_Karpas_genCR_Exp285_abs_LFC_d.sav', 'rb')))
Screen_name_list.append('Avana HCC1143 - Exp83')
Screen_data_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Avana_HCC1143_genCR_Exp227CEG3_abs_LFC_d.sav', 'rb')))
Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Avana_HCC1143_genCR_Exp227abs_LFC_d.sav', 'rb')))
Screen_name_list.append('Brunello A375')
Screen_data_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Brunello_CEG3_abs_LFC_d.sav', 'rb')))
Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Brunello_abs_LFC_d.sav', 'rb')))
Screen_name_list.append('IMP KBM7')
Screen_data_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/KBM7_medm1_abs_LFC_d.sav', 'rb')))
Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/KBM7_abs_LFC_d.sav', 'rb')))
Screen_name_list.append('IMP RKO')
Screen_data_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/RKO_medm1_abs_LFC_d.sav', 'rb')))
Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/RKO_abs_LFC_d.sav', 'rb')))
Screen_name_list.append('IMP MIApaca2')
Screen_data_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/MIAPACA2_medm1_abs_LFC_d.sav', 'rb')))
Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/MIAPACA2_abs_LFC_d.sav', 'rb')))

##################################################################################################################################
#################################################### MOUSE
'''
Species = 'Ms'
Screen_name_list = []
Screen_data_d_list = []
Screen_data_all_d_list = []
Screen_name_list.append('CrUMI mESC Data_2n')
Screen_data_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/CrUMI_2n_CEG3_abs_LFC_d.sav', 'rb')))
Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/CrUMI_2n_abs_LFC_d.sav', 'rb')))
Screen_name_list.append('CrUMI mESC Data_n')
Screen_data_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/CrUMI_n_CEG3_abs_LFC_d.sav', 'rb')))
Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/CrUMI_n_abs_LFC_d.sav', 'rb')))
Screen_name_list.append('IMP mESC Data_2n')
Screen_data_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Zub_2n_d18_CEG3_abs_LFC_d.sav', 'rb')))
Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Zub_2n_d18_abs_LFC_d.sav', 'rb')))
Screen_name_list.append('IMP mESC Data_n')
Screen_data_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Zub_n_d18_CEG3_abs_LFC_d.sav', 'rb')))
Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Zub_n_d18_abs_LFC_d.sav', 'rb')))

#print('load model_name')
#model_T1_AAw5_6prop_abs_LFC = pickle.load(open('/Users/georg.michlits/Desktop/Models/VBC_score/model_Nov_rel.sav','rb'))
# gn30 Ms is gn63 Hm
print('load Doench_with AApositional information')
D16_AA_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/D16_AA_d.sav', 'rb'))
Hart17_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Hart17_d.sav', 'rb'))
Scorer2_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Scorer2.0_d.sav', 'rb'))
D14_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/D14_d.sav', 'rb'))
print('load Doench_without AApositional information')
D16_woAA_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/D16_woAA_d.sav', 'rb'))
print('load inDelphi')
inDelphi = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/inDelphi_infr_d.sav', 'rb'))

inDelphi_fr_d = {}
for pos in inDelphi:
    inDelphi_fr_d[pos] = 1-inDelphi[pos]

#VBC_score_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/VBC_score_3IMPsc_CEG3_rel_v2_d.sav', 'rb'))
#VBC_score_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/VBC_score_KBM7_medm1_rel_v2_d.sav', 'rb'))
#VBC_score_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/VBC_score_KBM7_CEG3_rel_v2_d.sav', 'rb'))
#VBC_score_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/VBC_score_Nov_rel_v2_d.sav', 'rb'))
#VBC_score_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/VBC_score_3Nov_rel_v2_d.sav', 'rb'))
#VBC_score_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/VBC_score_3IMPsc_medm1_abs_v2_d.sav', 'rb'))
#VBC_score_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/VBC_score_3IMPsc_CEG3_abs_v2_d.sav', 'rb'))
#VBC_score_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/VBC_score_KBM7_medm1_abs_v2_d.sav', 'rb'))
#VBC_score_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/VBC_score_KBM7_CEG3_abs_v2_d.sav', 'rb'))

print('load other useful dictionaries')
Gene_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/' + Species + '/property_sav/Gene_d.sav', 'rb'))
CEG_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/' + Species + '/property_sav/CEG_pos_d.sav', 'rb'))
nonEG_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/' + Species + '/property_sav/nonEG_d.sav', 'rb'))
CN_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/' + Species + '/property_sav/CN_d.sav', 'rb'))



###########################################################################################
#############        Rerun the whole calculation for depleters only           #############
###########################################################################################

def det_dAUC_best_worst_depleters_only(Screen_data_d, VBC_score_d, D16_AA_d, data_name, model_name, print_plot=False):
    Structure_dict = {}
    Screen_data_d_cur = {}
    for pos in Screen_data_d:
        if pos in VBC_score_d:
            if pos in D16_AA_d:
                Screen_data_d_cur[pos] = Screen_data_d[pos]
                LFC = Screen_data_d[pos]
                gene = Gene_d[pos]
                VBC_score = VBC_score_d[pos]
                D16_score = D16_AA_d[pos]
                if not gene in Structure_dict:
                    Structure_dict[gene] = [[],[],[],[]]
                    ### [0] is pos, [1] doench score, [2] VBC score
                Structure_dict[gene][0].append(pos)
                Structure_dict[gene][1].append(D16_score)
                Structure_dict[gene][2].append(VBC_score)
                Structure_dict[gene][3].append(LFC)
    D16 = {}
    D16['best'] = {}
    D16['worst'] = {}
    VBC = {}
    VBC['best'] = {}
    VBC['worst'] = {}
    Rand = {}
    Rand['best'] = {}
    Rand['worst'] = {}
    Perfect = {}
    Perfect['best'] = {}
    Perfect['worst'] = {}
    n = 0
    guide_per_gene = []
    for gene in Structure_dict:
        n += 1
        triple_list = Structure_dict[gene]
        guide_per_gene.append(len(triple_list[0]))
        D16_best_pos = triple_list[0][triple_list[1].index(np.max(triple_list[1]))]
        Enrico_best_pos = triple_list[0][triple_list[2].index(np.max(triple_list[2]))]
        Perfect_best_pos = triple_list[0][triple_list[3].index(np.min(triple_list[3]))]
        D16_worst_pos = triple_list[0][triple_list[1].index(np.min(triple_list[1]))]
        Enrico_worst_pos = triple_list[0][triple_list[2].index(np.min(triple_list[2]))]
        Perfect_worst_pos = triple_list[0][triple_list[3].index(np.max(triple_list[3]))]
        Rand_worst_pos = triple_list[0][np.random.randint(len(triple_list[0]))]
        Rand_best_pos = triple_list[0][np.random.randint(len(triple_list[0]))]
        D16['best'][D16_best_pos] = 0
        D16['worst'][D16_worst_pos] = 0
        VBC['best'][Enrico_best_pos] = 0
        VBC['worst'][Enrico_worst_pos] = 0
        Rand['best'][Rand_best_pos] = 0
        Rand['worst'][Rand_worst_pos] = 0
        Perfect['best'][Perfect_best_pos] = 0
        Perfect['worst'][Perfect_worst_pos] = 0

    i = 0
    x = []
    EN_best_AC = []
    EN_best = 0
    EN_worst_AC = []
    EN_worst = 0
    D16_best_AC = []
    D16_best = 0
    D16_worst_AC = []
    D16_worst = 0
    Rand_best_AC = []
    Rand_best = 0
    Rand_worst_AC = []
    Rand_worst = 0
    Perfect_best_AC = []
    Perfect_best = 0
    Perfect_worst_AC = []
    Perfect_worst = 0
    for item in sorted(Screen_data_d.items(), key=lambda t: t[1]):
        i += 1
        x.append(i)
        pos = item[0]
        LFC = item[1]
        #print(pos + '\t' + str(LFC))
        if pos in VBC['best']:
            EN_best += 1/n
        if pos in VBC['worst']:
            EN_worst += 1/n
        EN_best_AC.append(EN_best)
        EN_worst_AC.append(EN_worst)
        if pos in D16['best']:
            D16_best += 1/n
        if pos in D16['worst']:
            D16_worst += 1/n
        D16_best_AC.append(D16_best)
        D16_worst_AC.append(D16_worst)
        if pos in Rand['best']:
            Rand_best += 1/n
        if pos in Rand['worst']:
            Rand_worst += 1/n
        Rand_best_AC.append(Rand_best)
        Rand_worst_AC.append(Rand_worst)
        if pos in Perfect['best']:
            Perfect_best += 1/n
        if pos in Perfect['worst']:
            Perfect_worst += 1/n
        Perfect_best_AC.append(Perfect_best)
        Perfect_worst_AC.append(Perfect_worst)

    if print_plot:
        plt.title('Train with' + model_name + 'test on' + data_name)
        plt.ylabel('AC_discovery')
        plt.xlabel('guide_ranking')
        plt.scatter(x, EN_best_AC, c='#FF00FF', s=1)
        plt.scatter(x, EN_worst_AC, c='#FF00FF', s=1)
        plt.scatter(x, D16_best_AC, c='#008000', s=1)
        plt.scatter(x, D16_worst_AC, c='#008000', s=1)
        plt.scatter(x, Rand_best_AC, c='black', s=1)
        plt.scatter(x, Rand_worst_AC, c='black', s=1)
        plt.scatter(x, Perfect_best_AC, c='grey', s=1)
        plt.scatter(x, Perfect_worst_AC, c='grey', s=1)
        plt.savefig(Analysis_name + '_' + data_name + '_' + model_name +'AUC.png')
        plt.close()
        plt.close()
        plt.close()
        plt.close()
        plt.close()
        plt.close()


    dAUC_VBC = (sum(EN_best_AC)-sum(EN_worst_AC))/i
    dAUC_D16 = (sum(D16_best_AC)-sum(D16_worst_AC))/i
    dAUC_rand = (sum(Rand_best_AC)-sum(Rand_worst_AC))/i
    dAUC_Perfect = (sum(Perfect_best_AC)-sum(Perfect_worst_AC))/i
    return (dAUC_VBC,dAUC_D16,dAUC_rand,dAUC_Perfect)


def det_dAUC_best_worst_all(Screen_data_all_d, VBC_score_d, D16_AA_d, data_name, model_name, print_plot=False):
    Structure_dict = {}
    Screen_data_d_cur = {}
    for pos in Screen_data_all_d:
        if pos in VBC_score_d:
            if pos in D16_AA_d:
                Screen_data_d_cur[pos] = Screen_data_d[pos]
                gene = Gene_d[pos]
                VBC_score = VBC_score_d[pos]
                D16_score = D16_AA_d[pos]
                if not gene in Structure_dict:
                    Structure_dict[gene] = [[],[],[]]
                    ### [0] is pos, [1] doench score, [2] VBC score
                Structure_dict[gene][0].append(pos)
                Structure_dict[gene][1].append(D16_score)
                Structure_dict[gene][2].append(VBC_score)
    D16 = {}
    D16['best'] = {}
    D16['worst'] = {}
    VBC = {}
    VBC['best'] = {}
    VBC['worst'] = {}
    Rand = {}
    Rand['best'] = {}
    Rand['worst'] = {}
    n = 0
    guide_per_gene = []
    for gene in Structure_dict:
        n += 1
        triple_list = Structure_dict[gene]
        guide_per_gene.append(len(triple_list[0]))
        D16_best_pos = triple_list[0][triple_list[1].index(np.max(triple_list[1]))]
        Enrico_best_pos = triple_list[0][triple_list[2].index(np.max(triple_list[2]))]
        D16_worst_pos = triple_list[0][triple_list[1].index(np.min(triple_list[1]))]
        Enrico_worst_pos = triple_list[0][triple_list[2].index(np.min(triple_list[2]))]
        Rand_worst_pos = triple_list[0][np.random.randint(len(triple_list[0]))]
        Rand_best_pos = triple_list[0][np.random.randint(len(triple_list[0]))]
        D16['best'][D16_best_pos] = 0
        D16['worst'][D16_worst_pos] = 0
        VBC['best'][Enrico_best_pos] = 0
        VBC['worst'][Enrico_worst_pos] = 0
        Rand['best'][Rand_best_pos] = 0
        Rand['worst'][Rand_worst_pos] = 0

    i = 0
    x = []
    EN_best_AC = []
    EN_best = 0
    EN_worst_AC = []
    EN_worst = 0
    D16_best_AC = []
    D16_best = 0
    D16_worst_AC = []
    D16_worst = 0
    Rand_best_AC = []
    Rand_best = 0
    Rand_worst_AC = []
    Rand_worst = 0
    for item in sorted(Screen_data_all_d.items(), key=lambda t: t[1]):
        i += 1
        x.append(i)
        pos = item[0]
        LFC = item[1]
        #print(pos + '\t' + str(LFC))
        if pos in VBC['best']:
            EN_best += 1/n
        if pos in VBC['worst']:
            EN_worst += 1/n
        EN_best_AC.append(EN_best)
        EN_worst_AC.append(EN_worst)
        if pos in D16['best']:
            D16_best += 1/n
        if pos in D16['worst']:
            D16_worst += 1/n
        D16_best_AC.append(D16_best)
        D16_worst_AC.append(D16_worst)
        if pos in Rand['best']:
            Rand_best += 1/n
        if pos in Rand['worst']:
            Rand_worst += 1/n
        Rand_best_AC.append(Rand_best)
        Rand_worst_AC.append(Rand_worst)


    if print_plot:
        plt.title('Train with' + model_name + 'test on' + data_name)
        plt.ylabel('AC_discovery')
        plt.xlabel('guide_ranking')
        plt.scatter(x, EN_best_AC, c='#FF00FF', s=1)
        plt.scatter(x, EN_worst_AC, c='#FF00FF', s=1)
        plt.scatter(x, D16_best_AC, c='#008000', s=1)
        plt.scatter(x, D16_worst_AC, c='#008000', s=1)
        plt.scatter(x, Rand_best_AC, c='black', s=1)
        plt.scatter(x, Rand_worst_AC, c='black', s=1)
        plt.savefig(data_name + '_' + model_name +'AUC.png')
        plt.close()
        plt.close()
        plt.close()
        plt.close()
        plt.close()
        plt.close()


    dAUC_VBC = (sum(EN_best_AC)-sum(EN_worst_AC))/i
    dAUC_D16 = (sum(D16_best_AC)-sum(D16_worst_AC))/i
    dAUC_rand = (sum(Rand_best_AC)-sum(Rand_worst_AC))/i
    return (dAUC_VBC,dAUC_D16,dAUC_rand)


def det_dAUC(Screen_data_all_d, CEG_d, nonEG_d):

    sum_total = 0
    sum_CEG = 0
    sum_nonEG = 0
    for pos_item in sorted(Screen_data_all_d.items(), key=lambda t: t[1]):
        pos = pos_item[0]
        LFC = pos_item[1]
        sum_total += 1
        if pos in CEG_d:
            sum_CEG += 1
        if pos in nonEG_d:
            sum_nonEG += 1
    AC_CEG = 0
    AC_nonEG = 0
    accum_CEG = 0
    accum_nonEG = 0
    x = []
    y_dia = []
    y_CEG = []
    y_nonEG = []
    i = 0
    for pos_item in sorted(Screen_data_all_d.items(), key=lambda t: t[1]):
        i += 1
        pos = pos_item[0]
        LFC = pos_item[1]
        x.append(i)
        y_dia.append(i/sum_total)
        if pos in CEG_d:
            AC_CEG += 1/sum_CEG
        if pos in nonEG_d:
            AC_nonEG += 1/sum_nonEG
        y_CEG.append(AC_CEG)
        accum_CEG += AC_CEG
        y_nonEG.append(AC_nonEG)
        accum_nonEG += AC_nonEG
    total_area = len(x)*1
    AUC_CEG = round(accum_CEG/total_area, 4)
    AUC_nonEG = round(accum_nonEG/total_area, 4)
    dAUC = round(AUC_CEG-AUC_nonEG, 4)
    return(dAUC,AUC_CEG,AUC_nonEG,sum_total,sum_CEG,sum_nonEG,x,y_dia,y_CEG,y_nonEG)


def det_CEG(Screen_data_all_d,CEG_d):
    CEG_list = []
    for pos in Screen_data_all_d:
        if pos in CEG_d:
            CEG_list.append(Screen_data_all_d[pos])
    CEG_depl = np.median(CEG_list)
    return(CEG_depl)


print('data\tmodel\tperfect\tVBC\tD16AA\trand\tBio\tinDeplphi\tD16woAA\trand\td_AUCess\tCEGdepl\tHart17sc\tD14score\tScorer2.0')

outfile = open(Analysis_name + '.txt','w')
outfile.write('data\tmodel\tperfect\tVBC\tD16AA\trand\tBio\tinDeplphi\tD16woAA\trand\td_AUCess\tCEGdepl\tHart17sc\tD14\tScorer2.0')

x = []
y_Bio = []
y_inD = []
y_D16 = []
y_D16wo = []
y_VBC = []
y_TR2 = []
y_rand = []
c_Bio = []
c_inD = []
c_D16 = []
c_D16wo = []
c_VBC = []
c_TR2 = []
c_rand = []
for i in range(len(Screen_name_list)):
    data_name = Screen_name_list[i]
    Screen_data_d = Screen_data_d_list[i]
    Screen_data_all_d = Screen_data_all_d_list[i]

    for model_name in model_name_list:

        #print('load VBC_score')
        VBC_score_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Apr1_VBC_score_'+model_name +'_d.sav', 'rb'))
        #print('load Bioscore')
        Bioscore_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Apr1_Bioscore_'+model_name +'_d.sav', 'rb'))
        dAUC_VBC,dAUC_D16,dAUC_rand1,dAUC_Perfect = det_dAUC_best_worst_depleters_only(Screen_data_d, VBC_score_d, D16_AA_d, data_name, model_name, print_plot=True)
        #print('VBC_score dAUC')
        #print(dAUC_VBC)
        #print('Doench2016_inclposition dAUC')
        #print(dAUC_D16)
        dAUC_Bio,dAUC_inDelphi,dAUC_rand2,dAUC_Perfect = det_dAUC_best_worst_depleters_only(Screen_data_d, Bioscore_d, inDelphi_fr_d, data_name, model_name)
        #print('Bioscore dAUC')
        #print(dAUC_Bio)
        #print('inDelphi_fr dAUC')
        #print(dAUC_inDelphi)
        dAUC_TR2,dAUC_D16woAA,dAUC_rand3,dAUC_Perfect = det_dAUC_best_worst_depleters_only(Screen_data_d, Hart17_d, D16_woAA_d, data_name, model_name)
        dAUC_D14,dAUC_Scorer2,dAUC_rand3,dAUC_Perfect = det_dAUC_best_worst_depleters_only(Screen_data_d, D14_d, Scorer2_d, data_name, model_name)
        #print('Doench2016_without_position dAUC')
        #print(dAUC_D16woAA)
        #print('mean of 3 random runs')
        dAUC_av_rand = np.mean([dAUC_rand1,dAUC_rand2,dAUC_rand3])
        #print(dAUC_av_rand)
        dAUC,AUC_CEG,AUC_nonEG,sum_total,sum_CEG,sum_nonEG,x,y_dia,y_CEG,y_nonEG = det_dAUC(Screen_data_all_d,CEG_d,nonEG_d)
        #print('essential_vs_nonessential d_AUC')
        #print(dAUC)
        CEG_depl = det_CEG(Screen_data_all_d, CEG_d)
        #print('CEG_median_depletion')
        #print(CEG_depl)

        print(data_name + ',' + model_name + ',' + str(dAUC_Perfect) + ',' + str(dAUC_VBC) + ',' + str(dAUC_D16) + ',' + str(dAUC_av_rand) + ',' + str(dAUC_Bio) + ',' + str(dAUC_inDelphi) + ',' + str(dAUC_D16woAA) + ',' +str(dAUC_av_rand) + ',' +str(dAUC) + ',' + str(CEG_depl)+ ',' + str(dAUC_TR2)+ ',' + str(dAUC_D14) + ',' + str(dAUC_Scorer2))
        line = (data_name + '\t' + model_name + '\t' + str(dAUC_Perfect) + '\t' + str(dAUC_VBC) + '\t' + str(dAUC_D16) + '\t' + str(dAUC_av_rand) + '\t' + str(dAUC_Bio) + '\t' + str(dAUC_inDelphi) + '\t' + str(dAUC_D16woAA) + '\t' +str(dAUC_av_rand) + '\t' +str(dAUC) + '\t' + str(CEG_depl) + '\t' + str(dAUC_TR2)+ '\t' + str(dAUC_D14) + '\t' + str(dAUC_Scorer2))
        outfile.write('\n' + line)
        x.append(CEG_depl)
        y_VBC.append(dAUC_VBC)
        c_VBC.append('c')
        y_Bio.append(dAUC_Bio)
        c_Bio.append('b')
        y_inD.append(dAUC_inDelphi)
        c_inD.append('orange')
        y_D16.append(dAUC_D16)
        c_D16.append('r')
        y_D16wo.append(dAUC_D16woAA)
        c_D16wo.append('m')
        y_TR2.append(dAUC_TR2)
        c_TR2.append('r')
        y_rand.append(dAUC_av_rand)
'''
plt.scatter(x,y_Bio,)
plt.scatter(x,y_TR2,)
plt.scatter(x,y_inD)
plt.scatter(x,y_rand)
plt.scatter(x,y_D16wo)
plt.savefig(Analysis_name +'all' '.pdf')
plt.close()
plt.close()
plt.close()
plt.close()
plt.close()
plt.scatter((x,y_VBC))
plt.scatter((x,y_D16))
plt.savefig(Analysis_name +'VBC_vs_ruleset2' '.pdf')
plt.close()
plt.close()
plt.scatter((x,y_VBC))
plt.scatter((x,y_TR2))
plt.savefig(Analysis_name +'VBC_vs_TR2' '.pdf')
'''