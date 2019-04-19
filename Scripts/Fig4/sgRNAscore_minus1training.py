__author__ = 'georg.michlits'

import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

#scaler = StandardScaler()
#scaler.fit(X_train)
#X_train = scaler.transform(X_train)
#X_test = scaler.transform(X_test)

Species = 'Hm'
Output_name = 'sgRNA_score_test_UMI'
Feat_plot_dir = Output_name + '_Feature_plots'
Coef_plot_dir = Output_name + '_Coef_plots'
Model_dir = Output_name + '_Model'
Fit_plot_dir = Output_name + '_Fit_plots'
os.mkdir(Feat_plot_dir)
outfile_mod3pfam = open(Feat_plot_dir + '/mod3pfam_data_out.txt', 'w')
os.mkdir(Coef_plot_dir)
os.mkdir(Model_dir)
os.mkdir(Fit_plot_dir)

#print('load phylo')
#phylo_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/phylo_v2_d.sav', 'rb'))
#print('load phast')
#phast_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/phast_v2_d.sav', 'rb'))
#print('load AA63')
#AA_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_cons_gn63_d.sav', 'rb'))
#print('load Doench')
#Hart17_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Hart17_d.sav', 'rb'))
#D16_woAA_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/D16_woAA_d.sav', 'rb'))
#tracrv2_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/tracrv2_mms_d.sav', 'rb'))
#tracrv2_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/tracrv2_mhm_d.sav', 'rb'))
#tracrv2_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/tracrv2_d.sav', 'rb'))
#IMPv2_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/IMP22_d.sav', 'rb'))
#print('load inDelphi')
#inDelphi_infr_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/inDelphi_infr_d.sav', 'rb'))
#print('load dist')
#distance_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/distance_d.sav', 'rb'))
#mod3_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/mod3_d.sav', 'rb'))
#print('load pfam')
#Pfam_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Pfam_all_d.sav', 'rb'))
#Pfam_domain_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Pfam_domain_d.sav', 'rb'))
#Pfam_family_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Pfam_family_d.sav', 'rb'))
#Pfam_repeat_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Pfam_repeat_d.sav', 'rb'))
#print('load T_streches')
#_4T_strech_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/_4T_strech_d.sav', 'rb'))
#over4T_strech_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/over4T_strech_d.sav', 'rb'))
#print('load cutting_frame')
#cut_frame_0_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/cut_frame_0_d.sav', 'rb'))
#cut_frame_1_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/cut_frame_1_d.sav', 'rb'))
#cut_frame_2_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/cut_frame_2_d.sav', 'rb'))
#print('load splice')
#Gene_direction_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Gene_direction_d.sav', 'rb'))
#exon_dist_start_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/exon_dist_start_d.sav', 'rb'))
#exon_dist_stop_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/exon_dist_stop_d.sav', 'rb'))

#print('load gaps')
#AA_gaps_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_gaps_d.sav', 'rb'))

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

CEG_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/' + Species + '/property_sav/CEG_pos_d.sav', 'rb'))
# VBC_hm
#data_abs_1 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/KBM7_medm1_rel_LFC_d.sav','rb'))
#data_abs_2 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/RKO_medm1_rel_LFC_d.sav','rb'))
#data_abs_3 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/MIAPACA2_medm1_rel_LFC_d.sav','rb'))

# Brunello
#data_abs_1 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Brunello_CEG3_rel_LFC_d.sav', 'rb'))

# TKOv3
#data_abs_5 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/TKOv3_HAP1d18_CEG3_rel_LFC_d.sav', 'rb'))

# Sanger
#Exp3 is HL60 Exp4 HT1080
#data_abs_1 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Yusa/Exp_data_thirdCEG/Exp3third_rel_LFC_d.sav', 'rb'))
#data_abs_2 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Yusa/Exp_data_thirdCEG/Exp4third_rel_LFC_d.sav', 'rb'))

# VBC Ms
#data_abs_1ms = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Zub_2n_d18_CEG3_rel_LFC_d.sav', 'rb'))
#data_abs_2ms = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Zub_n_d18_CEG3_rel_LFC_d.sav', 'rb'))

# UMI Ms
data_abs_1ms = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/CrUMI_2n_CEG3_rel_LFC_d.sav', 'rb'))
data_abs_2ms = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/CrUMI_n_CEG3_rel_LFC_d.sav', 'rb'))


i = 0
#for pos in data_abs:
#    print(pos + '\t' + str(data_abs[pos]))
#    i += 1
#    if i == 5:
#        break
for t in ['-Sa','-Br','-VB','-UM']:
    score_importance_train = {}
    score_importance_test = {}
    x_name = ['tracrv2','TKOv3','D16','comb','rand']
    Hm_screens = []#data_abs_1, data_abs_2, data_abs_3]#, data_abs_4, data_abs_5, data_abs_6, data_abs_7]
    Ms_screens = [data_abs_1ms, data_abs_2ms]# data_abs_3ms, data_abs_4ms]
    for x in x_name:
        score_importance_train[x] = []
        score_importance_test[x] = []
        i = 0
        i_training = 0
        #i_validation = 0
        guide_not_found = 0
        guide_found_wrong_AA_len = 0
        data = []
        c_list = []
        #data_valid = []
        target = []
        mod3_0 = []
        mod3_1 = []
        Pfam_0 = []
        Pfam_1 = []
        #target_valid = []
        Species = 'Hm'
        Hart17_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Hart17_d.sav', 'rb'))
        D16_woAA_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/D16_woAA_d.sav', 'rb'))
        if t == '-Sa':
            tracrv2_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/tracrv2_-Sanger_d.sav', 'rb'))
        if t == '-Br':
            tracrv2_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/tracrv2_-Brunello_d.sav', 'rb'))
        if t == '-VB':
            tracrv2_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/tracrv2_-VBC_d.sav', 'rb'))
        if t == '-UM':
            tracrv2_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/tracrv2_-UMI_d.sav', 'rb'))
        for data_set in Hm_screens:
            for pos in data_set:
                data_list = []
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
                i_training += 1
                #PHYLO
                ['tracrv2','D16','TKOv3','sgRNAscore','inDelphi','sgRNAscore+inDelphi','random']
                if x == 'tracrv2':
                    data_list.append(tracrv2_d[pos])
                if x == 'D16':
                    data_list.append(D16_woAA_d[pos])
                if x == 'TKOv3':
                    data_list.append(Hart17_d[pos])
                if x == 'comb':
                    #data_list.append(IMPv2_d[pos])
                    data_list.append(tracrv2_d[pos])
                    data_list.append(D16_woAA_d[pos])
                    data_list.append(Hart17_d[pos])
                if x == 'rand':
                    score = np.random.random()
                    data_list.append(score)

                #rand
                #data_list.append(np.random.random())
                # here you can insert if pos in CEG c.list append red
                c_list.append('k')
                #####################################################################################
                ## y - LFC data to fit
                data.append(data_list)
                target.append(-1*data_set[pos])
        Species = 'Ms'
        Hart17_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Hart17_d.sav', 'rb'))
        D16_woAA_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/D16_woAA_d.sav', 'rb'))
        if t == '-Sa':
            tracrv2_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/tracrv2_-Sanger_d.sav', 'rb'))
        if t == '-Br':
            tracrv2_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/tracrv2_-Brunello_d.sav', 'rb'))
        if t == '-VB':
            tracrv2_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/tracrv2_-VBC_d.sav', 'rb'))
        if t == '-UM':
            tracrv2_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/tracrv2_-UMI_d.sav', 'rb'))
        for data_set in Ms_screens:
            for pos in data_set:
                data_list = []
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
                i_training += 1
                #PHYLO
                ['tracrv2','D16','TKOv3','comb','rand']
                if x == 'tracrv2':
                    data_list.append(tracrv2_d[pos])
                if x == 'D16':
                    data_list.append(D16_woAA_d[pos])
                if x == 'TKOv3':
                    data_list.append(Hart17_d[pos])
                if x == 'comb':
                    #data_list.append(IMPv2_d[pos])
                    data_list.append(tracrv2_d[pos])
                    data_list.append(D16_woAA_d[pos])
                    data_list.append(Hart17_d[pos])
                if x == 'rand':
                    score = np.random.random()
                    data_list.append(score)
                #rand
                #data_list.append(np.random.random())
                # here you can insert if pos in CEG c.list append red
                c_list.append('k')
                #####################################################################################
                ## y - LFC data to fit
                data.append(data_list)
                target.append(-1*data_set[pos])
                i += 1

        print(str(i) + ' guides in dataset')
        #print(str(i_training) + ' guides in training dataset')
        #print(str(i_validation) + ' guides in validation dataset')
        #print(str(guide_found_wrong_AA_len) + ' guides found but on the edge 5pr or 3pr of CDS')
        print(str(guide_not_found) + ' guides skipped. missing transcript ID - no AA_cons scores')
        Feature_name_list = []
        ['blank','IMPv2','tracrv2','D16','TKOv3','sgRNAscore','random']
        if x == 'IMPv2':
            Feature_name_list.append('IMPv2')
        if x == 'tracrv2':
            Feature_name_list.append('tracrv2')
        if x == 'D16':
            Feature_name_list.append('Rule Set 2')
        if x == 'TKOv3':
            Feature_name_list.append('Hart2017')
        if x == 'comb':
            #Feature_name_list.append('IMPv2')
            Feature_name_list.append('tracrv2')
            Feature_name_list.append('Rule Set 2')
            Feature_name_list.append('Hart2017')
        if x == 'sgRNAscore-Zub':
            Feature_name_list.append('tracrv2')
            Feature_name_list.append('Rule Set 2')
            Feature_name_list.append('Hart2017')
        if x == 'inDelphi':
            Feature_name_list.append('inDelphi')
        if x == 'sgRNAscore+inDelphi':
            Feature_name_list.append('IMPv2')
            Feature_name_list.append('tracrv2')
            Feature_name_list.append('Rule Set 2')
            Feature_name_list.append('Hart2017')
            Feature_name_list.append('inDelphi')
        if x == 'rand':
            Feature_name_list.append('random')
        #Feature_name_list.append('random2')


        Target_name = 'log2 fold change'
        c_list_feat = ['k','k','k','k','k','k','k','k','k','k','k','k']

        data_ML = {'data': np.array(data),
                    'feature_names': np.array(Feature_name_list),
                    'target_names': np.array(Target_name),
                    'target': np.array(target)}

        print('Data input')
        print(data_ML['data'].shape)
        print('features')
        print(data_ML['feature_names'].shape)
        print('target')
        print(data_ML['target'].shape)

        i = 0
        for index, feature_name in enumerate(data_ML['feature_names']):
            plt.title(feature_name)
            plt.scatter(data_ML['target'], data_ML['data'][:,index], c=c_list, s=10, alpha=0.4)
            gradient, intercept, r_value, p_value, std_err = stats.linregress(data_ML['target'], data_ML['data'][:,index])
            plt.xlabel(Target_name, size=15)
            plt.ylabel(feature_name, size=15)
            mn = np.min(data_ML['target'])
            mx = np.max(data_ML['target'])
            x1 = np.linspace(mn, mx, 500)
            y1 = gradient*x1+intercept
            plt.plot(x1, y1, c_list_feat[i], linewidth=2)
            plt.plot(x1, y1, linewidth=2)
            i += 1
            plt.savefig(Feat_plot_dir + '/' + feature_name + '.png')
            plt.close()
            plt.close()
            plt.close()

        ###make mod3 plot
        outfile_mod3pfam.write('mod3_0\tmod3_1\tpfam_0\tpfam_1')
        n = max([len(mod3_0), len(mod3_1), len(Pfam_0), len(Pfam_1)])
        for i in range(n):
            if i < len(mod3_0):
                outfile_mod3pfam.write('\n' + str(mod3_0[i]))
            else:
                outfile_mod3pfam.write('\n')
            if i < len(mod3_1):
                outfile_mod3pfam.write('\t' + str(mod3_1[i]))
            else:
                outfile_mod3pfam.write('\t')
            if i < len(Pfam_0):
                outfile_mod3pfam.write('\t' + str(Pfam_0[i]))
            else:
                outfile_mod3pfam.write('\t')
            if i < len(Pfam_1):
                outfile_mod3pfam.write('\t' + str(Pfam_1[i]))
            else:
                outfile_mod3pfam.write('\t')

        RMS_list = []
        clf_score_list = []
        Test_score_list = []
        coef_list = []
        n_runs = 9
        print('current analysis: ' + x)
        for i in range(n_runs):
            X_train, X_test, y_train, y_test = train_test_split(data_ML['data'], data_ML['target'])

            #scaler = StandardScaler()
            #scaler.fit(X_train)

            #X_train = scaler.transform(X_train)
            #X_test = scaler.transform(X_test)

            clf = LinearRegression()
            clf.fit(X_train, y_train)

            print('clf_score train-test ' + str(clf.score(X_train,y_train)) + '-' + str(clf.score(X_test,y_test)))
            if not x == 'blank':
                print(x)
                score_importance_train[x].append(clf.score(X_train,y_train))
                score_importance_test[x].append(clf.score(X_test,y_test))
            else:
                score_importance_train[x].append(clf.score(X_train,y_train))
                score_importance_test[x].append(clf.score(X_test,y_test))

            clf_score_list.append('clf_score train-test ' + str(clf.score(X_train,y_train)) + '-' + str(clf.score(X_test,y_test)))

            predicted = clf.predict(X_test)
            expected = y_test
            plt.scatter(expected ,-1*predicted, alpha=0.4, s=10, c=c_list) # ULIS favourite style plt.scatter((np.log2(expected)), (1-(predicted)))
            plt.axis('tight')
            #plt.ylim(-4,-2)
            #plt.xlim(-10,0)
            #plt.plot([-10,0], [-10,0], '--k')
            plt.xlabel('True LFC-meassure', size=15)
            plt.ylabel('New score', size=15)
            gradient, intercept, r_value, p_value, std_err = stats.linregress(expected,-1*predicted)
            mn = np.min(expected)
            mx = np.max(expected)
            x1 = np.linspace(mn,mx,500)
            y1 = gradient*x1+intercept
            plt.plot(x1,y1, c='r', linewidth=2)
            plt.tight_layout()
            name = 'Fit_' + Output_name + '_' + str(i)
            plt.savefig(Fit_plot_dir + '/' + name + '.png')
            #plt.xlim(0,1)
            plt.close()

            print("RMS: %s" % np.sqrt(np.mean((predicted - expected) ** 2)))
            RMS_list.append(np.sqrt(np.mean((predicted - expected) ** 2)))

            Test_score_list.append(clf.score(X_test,y_test))
            coef_list.append(clf.coef_)
            plt.bar(Feature_name_list,clf.coef_)
            name = 'Coef_' + Output_name + '_' + str(i)
            plt.savefig(Coef_plot_dir + '/' + name + '.pdf')
            plt.close()
            model_name = 'model_' + Output_name + '_' + str(i)
            scaler_name = 'scaler_' + Output_name + '_' + str(i)
            pickle.dump(clf, open(Model_dir + '/' +model_name + '.sav', 'wb'))
            #pickle.dump(scaler, open(Model_dir + '/' +scaler_name + '.sav', 'wb'))

        i = Test_score_list.index(np.median(Test_score_list))
        print('Median performance model:')
        print('run_' + str(i))
        print(clf_score_list[i])
        print('load_model_' + str(i))
        median_model = pickle.load(open(Model_dir + '/model_' + Output_name + '_' + str(i) + '.sav', 'rb'))
        mm_predicted = median_model.predict(data_ML['data'])
        mm_expected = data_ML['target']
        median_model.score(data_ML['data'], data_ML['target'])


        plt.scatter(mm_expected, -1*mm_predicted, alpha=0.4, s=10, c=c_list) # ULIS favourite style plt.scatter((np.log2(expected)), (1-(predicted)))
        plt.axis('tight')
        #plt.ylim(-4,-2)
        #plt.xlim(-10,0)
        #plt.plot([-10,0], [-10,0], '--k')
        plt.xlabel('log2 fold change', size=15)
        plt.ylabel('New score', size=15)
        plt.tight_layout()
        name = 'MedianFit_' + Output_name + '_' + str(i)
        gradient, intercept, r_value, p_value, std_err = stats.linregress(mm_expected,-1*mm_predicted)
        mn = np.min(mm_expected)
        mx = np.max(mm_expected)
        x1 = np.linspace(mn,mx,500)
        y1 = gradient*x1+intercept
        plt.plot(x1,y1, c='r', linewidth=2)
        plt.savefig(name + '.png')
        #plt.xlim(0,1)
        plt.close()

        coef_reshape = []
        for element in Feature_name_list:
            coef_reshape.append([])

        for i in range(len(Feature_name_list)):
            for element in range(n_runs):
                coef_reshape[i].append(-1*coef_list[element][i])

        coef_mean_list = []
        coef_std_list = []
        for element in coef_reshape:
            coef_mean_list.append(np.mean(element))
            coef_std_list.append(np.std(element))
        plt.bar(Feature_name_list, coef_mean_list, color=c_list_feat, yerr=coef_std_list)
        name = 'Coef_mean_' + Output_name
        plt.savefig(name + '.png')
        plt.close()
        print(coef_mean_list)

        if x == 'blank':
            print('generate blank scores')
            blank_score_train = np.mean(score_importance_train['blank'])
            blank_score_test = np.mean(score_importance_test['blank'])
            print('blank_train: ' + str(blank_score_train) + '\t' + 'blank_test: ' + str(blank_score_test))


    x_name = []
    y_mean = []
    y_err = []
    ['IMPv2','tracrv2','D16','TKOv3','sgRNAscore','random']
    c_list = ['darkolivegreen','forestgreen','lawngreen','darkgreen','grey']
    #score_importance_train.pop('blank') #kick out the blank value (it is substracted from all other scores)
    print('printout_individual clf.scores on train_data_set')
    for x in score_importance_train:
        x_name.append(x)
        y_mean.append(np.mean(score_importance_train[x]))
        y_err.append(np.std(score_importance_train[x]))
        print(x + '\t' + str(y_mean) + '\t' + str(y_err))
    plt.bar(x_name, y_mean, color=c_list, yerr=y_err)
    name = 'Train_' + Output_name + t
    plt.tick_params(labelsize=16)
    #plt.ylim(0,0.12)
    plt.tight_layout()
    plt.savefig(name + '.pdf')
    plt.close()


    x_name = []
    y_mean = []
    y_err = []
    #score_importance_test.pop('blank') #kick out the blank value (it is substracted from all other scores)
    print('printout_individual clf.scores on test_data_set')
    for x in score_importance_test:
        x_name.append(x)
        y_mean.append(np.mean(score_importance_test[x]))
        y_err.append(np.std(score_importance_test[x]))
        print(x + '\t' + str(y_mean) + '\t' + str(y_err))

    plt.bar(x_name, y_mean, color=c_list, yerr=y_err)
    plt.tick_params(labelsize=16)
    #plt.ylim(0,0.12)
    plt.tight_layout()
    name = 'Test_' + Output_name + t
    plt.savefig(name + '.pdf')
    plt.close()
