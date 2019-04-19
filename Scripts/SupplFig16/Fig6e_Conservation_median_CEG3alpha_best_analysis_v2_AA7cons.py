__author__ = 'georg.michlits'

import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
from scipy import stats

Output_filename = 'Hm_KBM7_best_v2'
Species = 'Hm'

Gene_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Gene_d.sav', 'rb'))
Gene_conservation_gd = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA7Genes_gd.sav', 'rb'))
CN_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/CN_d.sav', 'rb'))
Bioscore_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Apr8_Bioscore_all_rel_v2_d.sav', 'rb'))

Screen_name_list = []
Screen_data_d_list = []
Screen_data_all_d_list = []
# Screen_name_list.append('CrUMI mESC Data_2n')
# Screen_data_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/CrUMI_2n_CEG3_abs_LFC_d.sav', 'rb')))
# Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/CrUMI_2n_abs_LFC_d.sav', 'rb')))
# Screen_name_list.append('CrUMI mESC Data_n')
# Screen_data_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/CrUMI_n_CEG3_abs_LFC_d.sav', 'rb')))
# Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/CrUMI_n_abs_LFC_d.sav', 'rb')))
# Screen_name_list.append('IMP mESC Data_2n')
# Screen_data_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Zub_2n_d18_CEG3_abs_LFC_d.sav', 'rb')))
# Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Zub_2n_d18_abs_LFC_d.sav', 'rb')))
# Screen_name_list.append('IMP mESC Data_n')
# Screen_data_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Zub_n_d18_CEG3_abs_LFC_d.sav', 'rb')))
# # Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Zub_n_d18_abs_LFC_d.sav', 'rb')))
Screen_name_list.append('IMP KBM7')
Screen_data_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/KBM7_medm1_abs_LFC_d.sav', 'rb')))
Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/KBM7_abs_LFC_d.sav', 'rb')))
# Screen_name_list.append('IMP RKO')
# Screen_data_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/RKO_medm1_abs_LFC_d.sav', 'rb')))
# Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/RKO_abs_LFC_d.sav', 'rb')))
# Screen_name_list.append('IMP MIApaca2')
# Screen_data_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/MIAPACA2_medm1_abs_LFC_d.sav', 'rb')))
# Screen_data_all_d_list.append(pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/MIAPACA2_abs_LFC_d.sav', 'rb')))

outfile = open(Output_filename + '.txt', 'w')

for screen_number in range(len(Screen_data_all_d_list)):
    screen_name = Screen_name_list[screen_number]
    screen = Screen_data_d_list[screen_number]
    i = 0
    k = 0
    l = 0
    Gene_dep_d = {}
    #Gene_dep_d[Gene][pos]= LFC
    for pos in screen:
        i += 1
        if pos in Gene_d:
            if CN_d[pos] == 1:
                gene = Gene_d[pos]
                if not gene in Gene_dep_d:
                    Gene_dep_d[gene] = {}
                    l += 1
                Gene_dep_d[gene][pos] = screen[pos]
                k += 1
    print(str(i) + 'guides in libary')
    print(str(k) + 'guides with genename and CN 1')
    print(str(l) + 'genes in analysis')

    # we now have:
    # Gene_dep_d[Gene][pos] = LFC
    # Gene_cons_d[Gene] = conservation_score
    # Determine average depletion per gene with at least 3 guides
    # Bin Genes in 3 groups of equal size based on conservation
    len_genes = len(Gene_conservation_gd)
    hi_med = int(len_genes/3)
    med_lo = int(len_genes/3*2)
    i = 0
    top_gene_depl = {}
    top_gene_depl['hi'] = {}
    top_gene_depl['me'] = {}
    top_gene_depl['lo'] = {}
    for item in sorted(Gene_conservation_gd.items(), key=lambda t: t[1], reverse=True):
        gene = item[0]
        conservation = item[1]
        #print(gene + '\t' + str(conservation))
        i += 1
        if i < hi_med:
            if gene in Gene_dep_d:
                if len(Gene_dep_d[gene]) >= 3:
                    LFC_list = []
                    pos_list = []
                    for pos in Gene_dep_d[gene]:
                        LFC_list.append(Gene_dep_d[gene][pos])
                        pos_list.append(pos)
                    top_gene_depl['hi'][gene] = [np.median(LFC_list),np.min(LFC_list),np.max(LFC_list),LFC_list,pos_list[LFC_list.index(np.min(LFC_list))],pos_list[LFC_list.index(np.max(LFC_list))],pos_list]
        elif i < med_lo:
            if gene in Gene_dep_d:
                if len(Gene_dep_d[gene]) >= 3:
                    LFC_list = []
                    pos_list = []
                    for pos in Gene_dep_d[gene]:
                        LFC_list.append(Gene_dep_d[gene][pos])
                        pos_list.append(pos)
                    top_gene_depl['me'][gene] = [np.median(LFC_list),np.min(LFC_list),np.max(LFC_list),LFC_list,pos_list[LFC_list.index(np.min(LFC_list))],pos_list[LFC_list.index(np.max(LFC_list))],pos_list]
        else:
            if gene in Gene_dep_d:
                if len(Gene_dep_d[gene]) >= 3:
                    LFC_list = []
                    pos_list = []
                    for pos in Gene_dep_d[gene]:
                        LFC_list.append(Gene_dep_d[gene][pos])
                        pos_list.append(pos)
                    print(LFC_list)
                    print(pos_list)
                    top_gene_depl['lo'][gene] = [np.median(LFC_list),np.min(LFC_list),np.max(LFC_list),LFC_list,pos_list[LFC_list.index(np.min(LFC_list))],pos_list[LFC_list.index(np.max(LFC_list))],pos_list]
    Hi_x = []
    Hi_y = []
    Hi_c = []
    Me_x = []
    Me_y = []
    Me_c = []
    Lo_x = []
    Lo_y = []
    Lo_c = []
    Hi_x1 = []
    Hi_y1 = []
    Hi_c1 = []
    Me_x1 = []
    Me_y1 = []
    Me_c1 = []
    Lo_x1 = []
    Lo_y1 = []
    Lo_c1 = []
    Hi_x2 = []
    Hi_y2 = []
    Hi_c2 = []
    Me_x2 = []
    Me_y2 = []
    Me_c2 = []
    Lo_x2 = []
    Lo_y2 = []
    Lo_c2 = []
    i = 0
    for item in sorted(top_gene_depl['hi'].items(), key=lambda t: t[1][1]):
        gene = item[0]
        gene_data = item[1]
        Hi_x.append(i)
        Hi_x1.append(i)
        Hi_x2.append(i)
        Hi_y.append(max([1-2**gene_data[0]]))
        Hi_y1.append(max([1-2**gene_data[2]]))
        Hi_y2.append(max([1-2**gene_data[1]]))
        outfile.write('\n' + 'hi' + '\t' + str(max([0,1-2**gene_data[2]])) + '\t' + str(max([0,1-2**gene_data[1]])))
        Hi_c2.append(Bioscore_d[gene_data[4]])
        Hi_c1.append(Bioscore_d[gene_data[5]])
        i+=1
    hi_genes = i
    print('high genes' + str(i))
    name = 'Highly Conserved genes'
    sctr = plt.scatter(Hi_x, Hi_y, c='k', s=3, linewidth=2)
    sctr = plt.scatter(Hi_x1, Hi_y1, c='r', s=1, alpha=0.5, linewidth=2)
    sctr = plt.scatter(Hi_x2, Hi_y2, c='b', s=1, alpha=0.5, linewidth=2)
    #plt.colorbar(sctr)
    #plt.title(name)
    #plt.xlabel('sgRNA depletion')
    #plt.ylabel('Genes sorted by depletion of stronges guide')
    #plt.savefig('high.pdf')
    for item in sorted(top_gene_depl['me'].items(), key=lambda t: t[1][1]):
        gene = item[0]
        gene_data = item[1]
        Me_x.append(i)
        Me_x1.append(i)
        Me_x2.append(i)
        Me_y.append(max([1-2**gene_data[0]]))
        Me_y1.append(max([1-2**gene_data[2]]))
        Me_y2.append(max([1-2**gene_data[1]]))
        outfile.write('\n' + 'me' + '\t' + str(max([0,1-2**gene_data[2]])) + '\t' + str(max([0,1-2**gene_data[1]])))
        Me_c2.append(Bioscore_d[gene_data[4]])
        Me_c1.append(Bioscore_d[gene_data[5]])
        i+=1
    me_genes = i-hi_genes
    print('medium genes' + str(me_genes))
    name = 'Medium Conserved genes'
    sctr = plt.scatter(Me_x, Me_y, c='k', s=3, linewidth=2)
    sctr = plt.scatter(Me_x1, Me_y1, c='r', s=1, alpha=0.5, linewidth=2)
    sctr = plt.scatter(Me_x2, Me_y2, c='b', s=1, alpha=0.5, linewidth=2)
    #plt.colorbar(sctr)
    #plt.title(name)
    #plt.xlabel('sgRNA depletion')
    #plt.ylabel('Genes sorted by depletion of stronges guide')
    #plt.savefig('medium.pdf')
    for item in sorted(top_gene_depl['lo'].items(), key=lambda t: t[1][1]):
        gene = item[0]
        gene_data = item[1]
        Lo_x.append(i)
        Lo_x1.append(i)
        Lo_x2.append(i)
        Lo_y.append(max([1-2**gene_data[0]]))
        Lo_y1.append(max([1-2**gene_data[2]]))
        Lo_y2.append(max([1-2**gene_data[1]]))
        outfile.write('\n' + 'lo' + '\t' + str(max([0,1-2**gene_data[2]])) + '\t' + str(max([0,1-2**gene_data[1]])))
        Lo_c2.append(Bioscore_d[gene_data[4]])
        Lo_c1.append(Bioscore_d[gene_data[5]])
        i+=1
    lo_genes = i-hi_genes-me_genes
    print('low genes' + str(lo_genes))
    name = 'High-medium-low conserved genes'
    sctr = plt.scatter(Lo_x, Lo_y, c='k', s=3, linewidth=2)
    sctr = plt.scatter(Lo_x1, Lo_y1, c='r', s=1, alpha=0.5, linewidth=2)
    sctr = plt.scatter(Lo_x2, Lo_y2, c='b', s=1, alpha=0.5, linewidth=2)
    #plt.colorbar(sctr)
    plt.title(name)
    plt.xlabel('3 Conservation Groups, Genes sorted by depletion of best sgRNA')
    plt.ylabel('best and worst sgRNA depletion %')
    plt.savefig(screen_name + 'low_linear_medv2_CEG2no0s.pdf')
    plt.close()

    pickle.dump(top_gene_depl,open('median_gene_depl_.sav','wb'))

'''
    #y = np.array(Hi_y)-np.array(Hi_y1)
    y = np.array(Hi_y1)
    x = np.array(Hi_y)
    sctr = plt.scatter(x,y, c='grey')
    gradient, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    mn = np.min(y)
    mx = np.max(y)
    x1 = np.linspace(mn, mx, 500)
    y1 = gradient*x1+intercept
    plt.title('High  Rscore ' + str(r_value) + '  p ' +str(p_value))
    plt.plot(x1, y1, c='r', linewidth=2)
    plt.xlabel('best sgRNA depletion')
    plt.ylabel('worst sgRNA depletion')
    plt.savefig('hi_linreg_v2.pdf')

    plt.close()
    # y = np.array(Me_y)-np.array(Me_y1)
    y = np.array(Me_y1)
    x = np.array(Me_y)
    sctr = plt.scatter(x,y, c='grey')
    gradient, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    mn = np.min(y)
    mx = np.max(y)
    x1 = np.linspace(mn, mx, 500)
    y1 = gradient*x1+intercept
    plt.title('Med   Rscore ' + str(r_value) + '  p ' +str(p_value))
    plt.xlabel('best sgRNA depletion')
    plt.ylabel('worst sgRNA depletion')
    plt.plot(x1, y1, c='r', linewidth=2)
    plt.savefig('me_linreg_v2.pdf')

    plt.close()
    #y = np.array(Lo_y)-np.array(Lo_y1)
    y = np.array(Lo_y1)
    x = np.array(Lo_y)
    sctr = plt.scatter(x,y, c='grey')
    gradient, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    mn = np.min(y)
    mx = np.max(y)
    x1 = np.linspace(mn, mx, 500)
    y1 = gradient*x1+intercept
    plt.title('Low  Rscore ' + str(r_value) + '  p ' +str(p_value))
    plt.xlabel('best sgRNA depletion')
    plt.ylabel('worst sgRNA depletion')
    plt.plot(x1, y1, c='r', linewidth=2)
    plt.savefig('lo_linreg_v2.pdf')
    plt.close()

    #only if best sgRNA about 80% depletion:


    #y = np.array(Hi_y)-np.array(Hi_y1)
    x = []
    y = []
    for i in range(len(Hi_y)):
        if Hi_y[i] > 0.8:
            x.append(Hi_y[i]*100)
            y.append(Hi_y1[i]*100)
    x = np.array(x)
    y = np.array(y)
    sctr = plt.scatter(x,y, c='grey')
    gradient, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    mn = np.min(y)
    mx = np.max(y)
    x1 = np.linspace(mn, mx, 500)
    y1 = gradient*x1+intercept
    plt.xlim(80,100)
    plt.ylim(0,100)
    plt.title('High  Rscore ' + str(r_value) + '  p ' +str(p_value))
    plt.plot(x1, y1, c='r', linewidth=2)
    plt.xlabel('best sgRNA depletion')
    plt.ylabel('worst sgRNA depletion')
    plt.savefig('hi_linreg_v3.pdf')
    plt.close()
    # y = np.array(Me_y)-np.array(Me_y1)
    x = []
    y = []
    for i in range(len(Me_y)):
        if Me_y[i] > 0.8:
            x.append(Me_y[i]*100)
            y.append(Me_y1[i]*100)
    x = np.array(x)
    y = np.array(y)
    sctr = plt.scatter(x,y, c='grey')
    gradient, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    mn = np.min(y)
    mx = np.max(y)
    x1 = np.linspace(mn, mx, 500)
    y1 = gradient*x1+intercept
    plt.xlim(80,100)
    plt.ylim(0,100)
    plt.title('Med   Rscore ' + str(r_value) + '  p ' +str(p_value))
    plt.xlabel('best sgRNA depletion')
    plt.ylabel('worst sgRNA depletion')
    plt.plot(x1, y1, c='r', linewidth=2)
    plt.savefig('me_linreg_v3.pdf')
    plt.close()
    #y = np.array(Lo_y)-np.array(Lo_y1)
    x = []
    y = []
    for i in range(len(Lo_y)):
        if Lo_y[i] > 0.8:
            x.append(Lo_y[i]*100)
            y.append(Lo_y1[i]*100)
    x = np.array(x)
    y = np.array(y)
    sctr = plt.scatter(x,y, c='grey')
    gradient, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    mn = np.min(y)
    mx = np.max(y)
    x1 = np.linspace(mn, mx, 500)
    y1 = gradient*x1+intercept
    plt.xlim(80,100)
    plt.ylim(0,100)
    plt.title('Low  Rscore ' + str(r_value) + '  p ' +str(p_value))
    plt.xlabel('best sgRNA depletion')
    plt.ylabel('worst sgRNA depletion')
    plt.plot(x1, y1, c='r', linewidth=2)
    plt.savefig('lo_linreg_v3.pdf')
'''
    #
    # window = 40
    # a = []
    # Hi_y1fl = []
    # for n in range(len(Hi_y1)):
    #     upper_pos_range = window/2
    #     lower_pos_range = window/2
    #     if n < window/2:
    #         upper_pos_range = 20 + (20-n)
    #         lower_pos_range = 20 - (20-n)
    #     if len(Hi_y1)-n < window/2:
    #         lower_pos_range = 20 + 20-len(Hi_y1)-n
    #         upper_pos_range = 20 + 20-len(Hi_y1)-n
    #     x = []
    #     for sgRNAdepl in










