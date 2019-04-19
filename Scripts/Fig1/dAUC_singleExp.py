__author__ = 'georg.michlits'

import pickle
import matplotlib.pyplot as plt

Species = 'Ms'
data_folder = ''    #enter URL of folder all .sav files will be loaded and analysed all plots stored in /plots data_name stored in 'NAME_AUC.txt'
data = 'Zub_ms_d6'
' + Species + '
Screen_data = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Zub_2n_d6_abs_LFC_d.sav','rb'))
CEG_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/' + Species + '/property_sav/CEG_pos_d.sav', 'rb'))
nonEG_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/' + Species + '/property_sav/nonEG_d.sav', 'rb'))

sum_total = 0
sum_CEG = 0
sum_nonEG = 0
for pos_item in sorted(Screen_data.items(), key=lambda t: t[1]):
    pos = pos_item[0]
    LFC = pos_item[1]
    sum_total += 1
    if pos in CEG_d:
        sum_CEG += 1
    if pos in nonEG_d:
        sum_nonEG += 1

print('total guides present: ' + str(sum_total))
print('total CEG guides present: ' + str(sum_CEG))
print('total non_Ess guides present: ' + str(sum_nonEG))

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

print('AUC_CEG: ' + str(AUC_CEG))
print('AUC_nonEG: ' + str(AUC_nonEG))
print('dAUC: ' + str(round(AUC_CEG - AUC_nonEG, 4)))

AUC_outfile = open(data + 'AUC.txt','w')

AUC_outfile.write('dAUC\tAUC_CEG\tAUC_nonEG\ttotal_guides\tCEG_guides\tnon_Ess_guides')
AUC_outfile.write('\n' + str(dAUC) + '\t' + str(AUC_CEG) + '\t' + str(AUC_nonEG) + '\t' + str(sum_total) + '\t' + str(sum_CEG) + '\t' + str(sum_nonEG))


plt.title('AUC ' + data)
plt.ylabel('AC_discovery')
plt.xlabel('guide_ranking')
plt.scatter(x, y_CEG, c='g', s=1)
plt.scatter(x, y_nonEG, c='r', s=1)
plt.scatter(x, y_dia, c='black', s=0.5)
plt.savefig('AUC ' + data + '.pdf')
plt.savefig('AUC ' + data + '.png')