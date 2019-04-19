__author__ = 'georg.michlits'

import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

Species = 'Ms'
Output_name = 'All_' + Species + 'VBC_guides'

#print('load Doench')
#tracrv2_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/tracrv2_d.sav', 'rb'))
D16_woAA_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/D16_woAA_d.sav', 'rb'))
#print('load inDelphi')
#inDelphi_infr_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/inDelphi_infr_d.sav', 'rb'))
#print('load Bioscore')
#Bioscore_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Apr1_Bioscore_3Nov_rel_v2_d.sav', 'rb'))

print('load genes')
Gene_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Gene_d.sav', 'rb'))
print('load VBC_sgRNA_score')
VBC_score_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Apr8_VBC_score_all_relraw_d.sav', 'rb'))
Bioscore_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Apr8_Bioscore_all_rel_v2_d.sav', 'rb'))
print('load distance')
distance_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/distance_d.sav', 'rb'))
print('load Doench')
D16_woAA_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/D16_woAA_d.sav', 'rb'))
tracrv2_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/tracrv2_d.sav', 'rb'))
print('load inDelphi')
inDelphi_infr_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/inDelphi_infr_d.sav', 'rb'))
print('load AA7aa')
AA_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA7_d.sav', 'rb'))
print('load copynumberfilter')
CN_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/CN_d.sav', 'rb'))
print('load AA_around cutsite')
AA_score_dm2 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_residue_d-2.sav', 'rb'))
AA_score_dm1 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_residue_d-1.sav', 'rb'))
AA_score_dv2 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_residue_dv2.sav', 'rb'))
AA_score_dp1 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_residue_d+1.sav', 'rb'))
AA_score_dp2 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_residue_d+2.sav', 'rb'))

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
'_' : 0}
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
nt30_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/' + Species + '/property_sav/nt30_d.sav', 'rb'))

Data_list_dict = {}
Structure_dict = {}
print('generatefeatures for model')
i = 0
guide_not_found = 0
data = []
for pos in D16_woAA_d:
    if pos in Gene_d:
        gene = Gene_d[pos]
        if gene not in Structure_dict:
            Structure_dict[gene] = {}
        Structure_dict[gene][pos] = VBC_score_d[pos]
    else:
        guide_not_found += 1
    i += 1
print(str(i) + ' guides in dataset')
print(str(guide_not_found) + ' guides skipped. missing transcript ID - no AA_cons scores')

tracr_v2_tf_d = {}
tr_list = []
for pos in tracrv2_d:
    tr_list.append(tracrv2_d[pos])
max_tr = max(tr_list)
min_tr = min(tr_list)
for pos in tracrv2_d:
    minmaxTrans = (tracrv2_d[pos]-min_tr)/(max_tr-min_tr)
    tracr_v2_tf_d[pos] = minmaxTrans

VBC_top50_file = open(Species + '_VBC_all.tsv','w')

VBC_top50_file.write('pos' + '\t' + 'gene' + '\t' + 'seq30N' + '\t' + 'VBC_score' + '\t' + 'Bioscore' + '\t' +  'sgRNA_activity_tracr_v2' + '\t' + 'sgRNA_activity_Doench2016_RuleSet2'
                     + '\t' + 'Frameshifts_inDelphi' + '\t' + 'perfect_cutsite_in_genome' +'\t' + 'distance_TSS=0_to_stop=1'  +'\t' 'targetes_AA_-2_to_+2' + '\t' + 'AA_identity_score' + '\t' + 'AA_cons7aa')
gene_count = 0
for gene in Structure_dict:
    gene_count += 1
    VBC_scores = {}
    for pos in Structure_dict[gene]:
        VBC_scores[pos] = VBC_score_d[pos]
    for item in sorted(VBC_scores.items(), key=lambda t: t[1], reverse=True):
        pos = item[0]
        VBC_top50_file.write('\n' + gene)
        VBC_top50_file.write('\t' + pos)
        VBC_top50_file.write('\t' + str(nt30_d[pos]))
        VBC_top50_file.write('\t' + str(round(VBC_score_d[pos],4)))
        VBC_top50_file.write('\t' + str(round(Bioscore_d[pos],4)))
        VBC_top50_file.write('\t' + str(round(tracr_v2_tf_d[pos],4)))
        VBC_top50_file.write('\t' + str(round(D16_woAA_d[pos],4)))
        if pos in inDelphi_infr_d:
            VBC_top50_file.write('\t' + str(round(1-inDelphi_infr_d[pos],4)))
        else:
            VBC_top50_file.write('\t' + str('NA'))
        VBC_top50_file.write('\t' + str(CN_d[pos]))
        VBC_top50_file.write('\t' + str(round(distance_d[pos],4)))
        if pos in AA_score_dv2:
            seq = ''
            seq += (AA_score_dm2[pos])
            seq += (AA_score_dm1[pos])
            seq += (AA_score_dv2[pos])
            seq += (AA_score_dp1[pos])
            seq += (AA_score_dp2[pos])
            data_list = []
            data_list.append(AA_score_m2[AA_score_dm2[pos]])
            data_list.append(AA_score_m1[AA_score_dm1[pos]])
            data_list.append(AA_score_0[AA_score_dv2[pos]])
            data_list.append(AA_score_p1[AA_score_dp1[pos]])
            data_list.append(AA_score_p2[AA_score_dp2[pos]])
            AA_score = np.mean(data_list)
            VBC_top50_file.write('\t' + str(seq))
            VBC_top50_file.write('\t' + str(round(AA_score,4)))
        else:
            VBC_top50_file.write('\t' + str('NA'))
            VBC_top50_file.write('\t' + str('NA'))
        if pos in AA_d:
            VBC_top50_file.write('\t' + str(round(AA_d[pos],4)))
        else:
            VBC_top50_file.write('\t' + str('NA'))