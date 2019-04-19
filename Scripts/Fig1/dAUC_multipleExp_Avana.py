__author__ = 'georg.michlits'

import pickle
import matplotlib.pyplot as plt
import os

Species = 'Hm'
# enter name that a new generated folder will have where data will be stored
new_data_folder_name = 'Avana_Hm'
os.mkdir(new_data_folder_name)
# enter URL of folder all .sav files will be loaded and analysed all plots stored in /plots data stored in 'NAME_AUC.txt'
data_folder = '/Users/georg.michlits/Desktop/Gen_data/Avana_cell_lines/all_LFCv1'

CEG_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/' + Species + '/property_sav/CEG_pos_d.sav', 'rb'))
nonEG_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/' + Species + '/property_sav/nonEG_d.sav', 'rb'))

def det_dAUC(Screen_data_d, CEG_d, nonEG_d):
    sum_total = 0
    sum_CEG = 0
    sum_nonEG = 0
    for pos_item in sorted(Screen_data_d.items(), key=lambda t: t[1]):
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
    for pos_item in sorted(Screen_data.items(), key=lambda t: t[1]):
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

AUC_outfile = open(new_data_folder_name + '_AUC.txt','w')
AUC_outfile.write('Exp_name\tdAUC\tAUC_CEG\tAUC_nonEG\ttotal_guides\tCEG_guides\tnon_Ess_guides')
for File_name in os.listdir(data_folder):
    if File_name.endswith('LFC_d.sav'):
        print('read_file ' + File_name)
        data_name = File_name
        Screen_data = pickle.load(open(data_folder + '/' + data_name,'rb'))
        dAUC,AUC_CEG,AUC_nonEG,sum_total,sum_CEG,sum_nonEG,x,y_dia,y_CEG,y_nonEG = det_dAUC(Screen_data, CEG_d, nonEG_d)
        print('read_file ' + File_name + '\t' + 'dAUC = ' + str(dAUC))
        AUC_outfile.write('\n' + data_name + '\t' + str(dAUC) + '\t' + str(AUC_CEG) + '\t' + str(AUC_nonEG) + '\t' + str(sum_total) + '\t' + str(sum_CEG) + '\t' + str(sum_nonEG))

        plt.title('AUC ' + data_name)
        plt.ylabel('AC_discovery')
        plt.xlabel('guide_ranking')
        plt.scatter(x, y_CEG, c='g', s=1)
        plt.scatter(x, y_nonEG, c='r', s=1)
        plt.scatter(x, y_dia, c='black', s=0.5)
        plt.savefig(new_data_folder_name + '/' + data_name + '.pdf')
        plt.savefig(new_data_folder_name + '/' + data_name + '.png')
        plt.close()