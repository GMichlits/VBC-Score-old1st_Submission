
__author__ = 'georg.michlits'

import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split

Species = 'Hm'
model_input = '3Nov_rel.sav'
Output_name = Species + '_make_Apr_Bio_score_' + model_input.split('.')[0] + ''
Feat_plot_dir = Output_name + '_Feature_plots'
Coef_plot_dir = Output_name + '_Coef_plots'
Model_dir = Output_name + '_Model'
Fit_plot_dir = Output_name + '_Fit_plots'
os.mkdir(Feat_plot_dir)
outfile_mod3pfam = open(Feat_plot_dir + '/mod3pfam_data_out.txt','w')
os.mkdir(Coef_plot_dir)
os.mkdir(Model_dir)
os.mkdir(Fit_plot_dir)

print('load model')
model_Bioscore = pickle.load(open('/Users/georg.michlits/Desktop/Models/Bioscore/Apr1model_' + model_input,'rb'))
print('load phylo')
phylo_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/phylo21_d.sav', 'rb'))
print('load phast')
phast_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/phast21_d.sav', 'rb'))
print('load AA63')
if Species == 'Ms':
    AA_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA7_d.sav', 'rb'))
else:
    AA_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA7_d.sav', 'rb'))
print('load Doench')
D16_woAA_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/tracrv2_d.sav', 'rb'))
print('load inDelphi')
inDelphi_infr_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/inDelphi_infr_d.sav', 'rb'))
print('load dist')
distance_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/distance_d.sav', 'rb'))
mod3_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/mod3_d.sav', 'rb'))
print('load pfam')
#Pfam_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Pfam_all_d.sav', 'rb'))
Pfam_domain_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Pfam_domain_d.sav', 'rb'))
Pfam_family_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Pfam_family_d.sav', 'rb'))
#Pfam_repeat_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Pfam_repeat_d.sav', 'rb'))
print('load copynumberfilter')
CN_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/CN_d.sav', 'rb'))
#print('load T_streches')
#_4T_strech_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/_4T_strech_d.sav', 'rb'))
#over4T_strech_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/over4T_strech_d.sav', 'rb'))
#print('load cutting_frame')
#cut_frame_0_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/cut_frame_0_d.sav', 'rb'))
#cut_frame_1_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/cut_frame_1_d.sav', 'rb'))
#cut_frame_2_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/cut_frame_2_d.sav', 'rb'))
print('load splice')
Gene_direction_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Gene_direction_d.sav', 'rb'))
exon_dist_start_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/exon_dist_start_d.sav', 'rb'))
exon_dist_stop_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/exon_dist_stop_d.sav', 'rb'))
print('load AAtype')
AA_score_dm2 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_residue_d-2.sav', 'rb'))
AA_score_dm1 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_residue_d-1.sav', 'rb'))
AA_score_dv2 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_residue_dv2.sav', 'rb'))
AA_score_dp1 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_residue_d+1.sav', 'rb'))
AA_score_dp2 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_residue_d+2.sav', 'rb'))
AA_type_cut = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_type_cut_d.sav', 'rb'))
AA_type_down = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_type_down_d.sav', 'rb'))
AA_type_up = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_type_up_d.sav', 'rb'))

###################
# NOTES
#   All dictionaries with _d contain "NAME_Property"_d = {pos : score}
#   All scores are between 0 and 1
#   pos is a guide ID string: chr12:85822107-85822137(+)
#   Dictionaries AA_cons_gn30_d has option of choosing window size "NAME_Property"_d = {pos : {1: score,
#                                                                                         3: score,
#                                                                                         5: score,.....}
#
###################
'''
Gene_direction_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/Gene_direction_d.sav', 'rb'))
exon_dist_start_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/exon_dist_start_d.sav', 'rb'))
exon_dist_stop_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/exon_dist_stop_d.sav', 'rb'))
Gene_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/Gene_d.sav', 'rb'))
UMI_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/UMI_d.sav', 'rb'))
distance_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/distance_d.sav', 'rb'))
mod3_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/mod3_d.sav', 'rb'))
Pfam_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/Pfam_all_d.sav', 'rb'))
Pfam_domain_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/Pfam_domain_d.sav', 'rb'))
Pfam_family_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/Pfam_family_d.sav', 'rb'))
Pfam_repeat_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/Pfam_repeat_d.sav', 'rb'))
ATG_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/ATG_d.sav', 'rb'))
_4T_strech_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/_4T_strech_d.sav', 'rb'))
over4T_strech_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/over4T_strech_d.sav', 'rb'))
transcript_ID_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/transcriptID_d.sav', 'rb'))
cut_frame_0_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/cut_frame_0_d.sav', 'rb'))
cut_frame_1_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/cut_frame_1_d.sav', 'rb'))
cut_frame_2_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/cut_frame_2_d.sav', 'rb'))
'''

#data_abs = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Brunello/Br_medm0.5_abs_LFC_d.sav','rb'))


all_AA = 'RHKDESTNQGAVILMFYWCP'
AA_score_m2 = {
'H':2.057008194,
'R':1.666008011,
'K':-0.631652433,
'D':-3.913610482,
'E':0.266781658,
'S':-3.152017762,
'T':2.505081059,
'N':-3.246495686,
'Q':-2.353445559,
'G':1.135708375,
'A':2.691284761,
'V':2.87413591,
'I':1.061684877,
'L':-0.604897702,
'M':-5.318803846,
'F':-2.214438057,
'Y':3.896367429,
'W':-0.730613099,
'C':8.886985888,
'P':-4.584034852,
'_':0}

AA_score_m1 = {
'H':0.248589167,
'R':1.608083504,
'K':-5.112688192,
'D':-3.713168458,
'E':1.277203186,
'S':-3.212790961,
'T':1.42280977,
'N':-2.495366944,
'Q':-2.33151425,
'G':-4.008199512,
'A':-4.316121484,
'V':-2.147544408,
'I':5.065845529,
'L':3.407003179,
'M':2.453279017,
'F':3.711715204,
'Y':12.00970109,
'W':-0.513787269,
'C':0.539820824,
'P':-1.817886089,
'_':0}

AA_score_p2 = {
'H':0.32397073,
'R':4.08757342,
'K':-8.369304206,
'D':-1.8599482,
'E':-2.663715943,
'S':-0.569723883,
'T':-1.699648991,
'N':-4.654929894,
'Q':-1.956412476,
'G':-2.005119053,
'A':2.673550002,
'V':3.531719149,
'I':1.591998958,
'L':3.869458891,
'M':-3.455773192,
'F':7.710046533,
'Y':6.560291786,
'W':5.741555436,
'C':3.800124529,
'P':-5.174282001,
'_':0}


AA_score_p1 = {
'H':-0.430372016,
'R':1.227731271,
'K':-4.162989207,
'D':-1.958772238,
'E':-1.255202252,
'S':-4.037535968,
'T':-7.086967292,
'N':-8.33702616,
'Q':-4.272260102,
'G':-0.444477914,
'A':-3.416425835,
'V':5.243732546,
'I':7.449137764,
'L':4.031869555,
'M':2.290774347,
'F':8.006760249,
'Y':8.409137581,
'W':5.54935109,
'C':9.564710848,
'P':-5.184584537,
'_':0}


AA_score_0 = {
'H' : 3.40143101,
'R' : 1.08052375,
'K' : -4.453883147,
'D' : -7.705176369,
'E' : -2.049637515,
'S' : -5.176921052,
'T' : 3.145321586,
'N' : -5.88211047,
'Q' : -2.337120609,
'G' : 3.364400709,
'A' : -5.766952958,
'V' : 3.855088089,
'I' : 3.994762418,
'L' : 2.750942054,
'M' : 1.06067543,
'F' : -1.795480355,
'Y' : 8.083789415,
'W' : 13.8994159,
'C' : 6.413725606,
'P' : -1.248312289,
'_':0}

phylo_list = []
phast_list = []
AA_list = []
Pfam_family_blank = 0.5
Pfam_domain_blank = 0.5
distance_blank = 0.5
mod3_blank = 0.5
exlen_list = []
SA_blank = 21
SD_blank = 22
AA_score = 0
AA_letters = 'RHKDESTNQGAVILMFYWCP'
AA_letters = '00000000000000000000'

for pos in phylo_d:
    phylo_list.append(phylo_d[pos])
median_phylo_list = np.median(phylo_list)
print('median phylo ' + str(median_phylo_list))
for pos in phast_d:
    phast_list.append(phast_d[pos])
median_phast_list = np.median(phast_list)
print('median phast ' + str(median_phast_list))
for pos in AA_d:
    AA_list.append(AA_d[pos])
median_AA_list = np.median(AA_list)
print('median phylo ' + str(median_AA_list))
for pos in exon_dist_stop_d:
    exlen_list.append(exon_dist_stop_d+exon_dist_start_d)
median_exlen_list = np.median(exlen_list)
print('median exlen ' + str(median_exlen_list))



VBC_crude_score_d = {}
Bioscore_crude_score_d = {}
VBC_list = []
Bioscore_list = []
total_guides = 0
guide_excluded = 0
for pos in D16_woAA_d:
    if pos in AA_d and pos in inDelphi_infr_d and pos in Gene_direction_d and pos in AA_score and pos in exon_dist_start_d:
        data_list = []
        total_guides += 1
        #if np.random.random() <= 0.1:
        #    i_validation += 1
        #    U6_list = AA_U6PAM_d[pos]['U6']
        #    PAM_list = AA_U6PAM_d[pos]['PAM']
        #    if len(U6_list)+len(PAM_list) == 21:
        #        U6_rev = []
        #        for x in reversed(U6_list):
        #            U6_rev.append(x)
        #        list = U6_rev + PAM_list[1:]
        #        data_valid.append(list)
        #        target_valid.append(2**target_data_file[pos])
        #    else:
        #        guide_found_wrong_AA_len += 1
        #else:
        #####################################################################################
        ##                    x Features  for optimization
        phylo = phylo_d[pos]
        data_list.append(phylo)
        score = phast_d[pos]
        data_list.append(score)
        score = Pfam_domain_d[pos]
        data_list.append(score)
        score = Pfam_family_d[pos]
        data_list.append(score)
        score = distance_d[pos]
        data_list.append(score)
        data_list.append(mod3_d[pos])
        data_list.append((exon_dist_start_d[pos]+exon_dist_stop_d[pos]))
        if Gene_direction_d[pos] == '+':
            score = exon_dist_start_d[pos]
            if score > 20:
                data_list.append(21)
            else:
                data_list.append(score)
        else:
            score = exon_dist_stop_d[pos]
            if score > 20:
                data_list.append(21)
            else:
                data_list.append(score)
        if Gene_direction_d[pos] == '-':
            score = exon_dist_start_d[pos]
            if score > 20:
                data_list.append(21)
            else:
                data_list.append(score)
        else:
            score = exon_dist_stop_d[pos]
            if score > 20:
                data_list.append(21)
            else:
                data_list.append(score)
        score = AA_d[pos]
        data_list.append(score)
        data_list.append(AA_score_m2[AA_score_dm2[pos]])
        data_list.append(AA_score_m1[AA_score_dm1[pos]])
        data_list.append(AA_score_0[AA_score_dv2[pos]])
        data_list.append(AA_score_p1[AA_score_dp1[pos]])
        data_list.append(AA_score_p2[AA_score_dp2[pos]])
        for letter in all_AA:
            UP_down_range = 6 # 6 corres[onds to window 13: cutsite 1 + 6 up  +6 down.
            if letter in AA_type_cut[pos]:
                data_list.append(1)
            elif letter in AA_type_up[pos][0:UP_down_range]:
                data_list.append(1)
            elif letter in AA_type_down[pos][0:UP_down_range]:
                data_list.append(1)
            else:
                data_list.append(0)
        test = 'ok'
        for element in data_list:
            if np.isnan(element):
                test = 'nan'
        if not 1 >= AA_d[pos] >= 0:
            test = 'nan'
        if test == 'ok':
            Bioscore = round(model_Bioscore.predict([data_list])[0],5)*-1
            Bioscore_crude_score_d[pos] = Bioscore
            Bioscore_list.append(Bioscore)
        #score = tracrv2_d[pos]
        #data_list.append(score)
        #score = 1-inDelphi_infr_d[pos]
        #data_list.append(score)
        #test = 'ok'
        #for element in data_list:
        #    if np.isnan(element):
        #        test = 'nan'
        #if not 1 >= AA_d[pos][7] >= 0:
        #    test = 'nan'
        #if test == 'ok':
        #    VBC_score = round(model_Bioscore.predict([data_list])[0],5)*-1
        #    VBC_crude_score_d[pos] = VBC_score
        #    VBC_list.append(VBC_score)
        #print('VBC ' + str(VBC_score) + '\t' + 'Bioscore ' + str(Bioscore))
    if total_guides % 100000 == 0:
        print(str(total_guides/1000000) + ' million guides processed')

print('total guides: ' +str(total_guides))
print('guides_excluded: ' +str(guide_excluded))

#VBC_max = np.max(VBC_list)
#VBC_min = np.min(VBC_list)
Biosc_max = np.max(Bioscore_list)
Biosc_min = np.min(Bioscore_list)

#plt.scatter(VBC_list,Bioscore_list,alpha=0.2)
#plt.xlabel = 'VBC_score_crude'
#plt.ylabel = 'Bio_score_crude'
#plt.savefig(Output_name + '_VBCscBio_crudeScores.png')

#print(VBC_max)
#print(VBC_min)
print(Biosc_max)
print(Biosc_min)

#VBC_score_d = {}
#for pos in VBC_crude_score_d:
#    VBC_scaled = (VBC_crude_score_d[pos]-VBC_min)/(VBC_max-VBC_min)
#    if VBC_scaled > 0.6:
#        VBC_scaledv2 = (VBC_scaled-0.6)/0.4*0.9+0.1
#    else:
#        VBC_scaledv2 = (VBC_scaled)/0.6*0.1
#    VBC_score_d[pos] = VBC_scaledv2

Bioscore_d = {}
for pos in Bioscore_crude_score_d:
    Bioscore_scaled = (Bioscore_crude_score_d[pos]-Biosc_min)/(Biosc_max-Biosc_min)
    if Bioscore_scaled > 0.6:
        Bio_scaledv2 = (Bioscore_scaled-0.6)/0.4*0.9+0.1
    else:
        Bio_scaledv2 = (Bioscore_scaled)/0.6*0.1
    Bioscore_d[pos] = Bio_scaledv2

#pickle.dump(VBC_score_d, open('HM_XXX_VBC_score_' + model_input.split('.')[0] +'_v2_d.sav','wb'))
pickle.dump(Bioscore_d, open(Species + '_Apr1_Bioscore_' + model_input.split('.')[0] + '_v2_d.sav','wb'))