__author__ = 'georg.michlits'


import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
from scipy import stats

def revcomp(dna, reverse=True, complement=True):
    DNA = dna.upper()
    a = 'ACGTRYKMSWBDHVN'
    b = 'TGCAYRMKSWVHDBN'
    dict = {a[i]:b[i] for i in range (15)}
    if reverse:
        DNA = reversed(DNA)
    if complement:
        result = [dict[i] for i in DNA]
    else: result = [i for i in DNA]
    return ''.join(result)

Species = 'Hm'
Outfile = open('top6_selected_hm_VBC_sgRNAs.tsv','w')

print('load VBC_sgRNA_score, etc.')
VBC_score_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Apr8_VBC_score_all_relraw_d.sav', 'rb'))
Bioscore_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Apr8_Bioscore_all_rel_v2_d.sav', 'rb'))
D16_woAA_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+ Species +'/property_sav/D16_woAA_d.sav', 'rb'))
tracrv2_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+ Species +'/property_sav/tracrv2_d.sav', 'rb'))
inDelphi_infr_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+ Species +'/property_sav/inDelphi_infr_d.sav', 'rb'))

AA_score_dm2 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_residue_d-2.sav', 'rb'))
AA_score_dm1 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_residue_d-1.sav', 'rb'))
AA_score_dv2 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_residue_dv2.sav', 'rb'))
AA_score_dp1 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_residue_d+1.sav', 'rb'))
AA_score_dp2 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_residue_d+2.sav', 'rb'))

print('load TTTTT-strech and CN dictinaries')
over4T_strech_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/over4T_strech_d.sav', 'rb'))
CN_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/CN_d.sav', 'rb'))
nt30_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/' + Species + '/property_sav/nt30_d.sav', 'rb'))
print('load Gene names')
Gene_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/' + Species + '/property_sav/Gene_d.sav', 'rb'))
print('load distance')
distance_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/' + Species + '/property_sav/distance_d.sav', 'rb'))

top50_OT_file = open(Species+'/offtarget.txt','r')

print('read in sgRNAs and penalize for copynumber and offtargets, etc.')
sgRNAs = {}
natural_g_dict = {18:{},19:{},20:{}}
OT90_d = {}
# sgRNAs[gene][pos] = VBC_picking_score
for line in top50_OT_file:
    col = line.rstrip('\n').split('\t')
    pos = col[0]
    OT90 = int(col[2])
    OT90_d[pos] = OT90
    CN = CN_d[pos]
    VBC = VBC_score_d[pos]
    Ts = over4T_strech_d[pos]
    seq = nt30_d[pos].upper()
    # determine if it is possible to make a guide that starts with natural G (for length 20bp 19 bp or 18bp):
    if seq[4] == 'G':
        naturalG_20 = 1
        natural_g_dict[20][pos]=0
    else:
        naturalG_20 = 0
    if seq[5] == 'G':
        naturalG_19 = 1
        natural_g_dict[19][pos]=0
    else:
        naturalG_19 = 0
    if seq[6] == 'G':
        naturalG_18 = 1
        natural_g_dict[18][pos]=0
    else:
        naturalG_18 = 0
    VBC_score_picking = VBC
    # penality for copy number (number of perfect target sites in the genome):
    if CN == 2:
        VBC_score_picking = VBC_score_picking * 0.70
    elif CN >= 3:
        VBC_score_picking = VBC_score_picking * 0.10
    # penality for copy number (OT90 is an estimate of likely off-target sites in the genome):
    if OT90 == 1:
        VBC_score_picking = VBC_score_picking * 0.80
    elif OT90 == 2:
        VBC_score_picking = VBC_score_picking * 0.60
    elif OT90 == 3:
        VBC_score_picking = VBC_score_picking * 0.40
    elif OT90 == 4:
        VBC_score_picking = VBC_score_picking * 0.20
    elif 4 < OT90 < 10:
        VBC_score_picking = VBC_score_picking * 0.10
    elif OT90 >= 10:
        VBC_score_picking = VBC_score_picking * 0.01
    # penality for TTTTT strech (5T in a row act as PolIII terminator sequence and should be avoided):
    if Ts == 1:
        VBC_score_picking = VBC_score_picking * 0.01
    if naturalG_19 == 1:
        VBC_score_picking = VBC_score_picking *1.2
    elif naturalG_18 == 1:
        VBC_score_picking = VBC_score_picking *1.1
    elif naturalG_20 == 1:
        VBC_score_picking = VBC_score_picking *1.0
    else:
        VBC_score_picking = VBC_score_picking *0.66
    genename = Gene_d[pos]
    if not genename in sgRNAs:
        sgRNAs[genename] = {}
    sgRNAs[genename][pos] = VBC_score_picking
    #print(genename + '\t' + str(VBC) + '\t' +str(VBC_score_picking) + '\t' + str(OT90) + '\t' + str(CN) + '\t' + str(Ts))


Outfile.write('genename' + '\t' + 'chr_position' + '\t' + '30bp_sequence'+ '\t' + 'VBC_score_0=bad_1=good' + '\t' + 'cutting_AA_seq_-2_to+2' + '\t' +
                          'distance_TSS=0_to_stop=1' + '\t' + 'oligotype' + '\t' + 'FW_oligo_name' + '\t' + 'FW_oligo_seq' + '\t' +
                          'RV_oligo_name' + '\t' + 'RV_oligo_seq' + '\t' + 'Bioscore_0=bad_1=good' + '\t' +
                          'sgRNA_activity_tracrv2_-0.5=bad_0.5=good' + '\t' + 'sgRNA_activity_Doench_0=bad_1=good' + '\t' + 'inDelphi_fraction_predicted_frameshifts' + '\t' +
                          'copy_number=perfect_match_target_sites_in_entire_genome' + '\t' + 'Off-targets_number_of_likely_OTs_in_codingDNAsequence')

print('run selection process pick 6 sgRNA per gene')
# here we run the selection for every gene 6 times (to get 6 guides) one by one we modify the vbc_score_picking values
# depending on proximity to regions that have been target with previously every selected guide varies the vbc_score_picking_values.
for gene in sorted(sgRNAs):
    #hitzones are entered as lists [start,stop]
    #print(gene)
    hitzone = []
    for i in range(6):
        for item in sorted(sgRNAs[gene].items(), key = lambda t: t[1], reverse=True):
            pos = item[0]
            VBC_score_picking = item[1]
            dist = distance_d[pos]
            hit_type = '1st'
            # check if it is a 'second hit within a 10% range of a previous hit zone in the protein
            for element in hitzone:
                if element[0] < dist < element[1]:
                    hit_type = '2nd'
                    # extend the size of the zone depending on the new hit
                    new_zone = [min(element[0],dist-0.05),max(element[1],dist+0.05)]
                    hitzone[hitzone.index(element)] = new_zone
                    break
            # remove selected sgRNA from dictionary
            sgRNAs[gene].pop(pos)

            ####################################################
            ##########     get guide properties and write result
            ####################################################
            #print(pos + '\t' + str(dist))
            if pos in natural_g_dict[19]:
                guide_sequence = nt30_d[pos][5:24]
                oligotype = '19nt_naturalG'
            elif pos in natural_g_dict[20]:
                guide_sequence = nt30_d[pos][4:24]
                oligotype = '20nt_naturalG'
            elif pos in natural_g_dict[18]:
                guide_sequence = nt30_d[pos][6:24]
                oligotype = '18nt_naturalG'
            else:
                guide_sequence = 'G' + nt30_d[pos][4:24]
                oligotype = '20nt_extra'

            FW_oligo = 'CACC' + guide_sequence
            RV_oligo = revcomp(guide_sequence + 'GTTT')

            if pos in AA_score_dv2:
                AAseq = ''
                AAseq += (AA_score_dm2[pos])
                AAseq += (AA_score_dm1[pos])
                AAseq += (AA_score_dv2[pos])
                AAseq += (AA_score_dp1[pos])
                AAseq += (AA_score_dp2[pos])
            else:
                AAseq = 'N.A.'

            if pos in inDelphi_infr_d:
                frameshifts = str(1-inDelphi_infr_d[pos])
            else:
                frameshifts = 'N.A.'

            Outfile.write('\n' + gene + '\t' + pos + '\t' + nt30_d[pos] + '\t' + str(VBC_score_d[pos]) + '\t' + AAseq + '\t' +
                          str(distance_d[pos]) + '\t' + oligotype + '\t' + gene +'_sg'+str(i+1)+'_FW' + '\t' + FW_oligo + '\t' +
                          gene +'_sg'+str(i+1)+'_RV' + '\t' + RV_oligo + '\t' + str(Bioscore_d[pos]) + '\t' +
                          str(tracrv2_d[pos]) + '\t' + str(D16_woAA_d[pos]) + '\t' + frameshifts + '\t' +
                          str(CN_d[pos]) + '\t' + str(OT90_d[pos]))

            ####################################################
            ##########     continue with guide selection
            ####################################################

            # rescore all other guides based on position relative to new hitzone (distinguish 2nd or 1st hit)
            if hit_type == '2nd':
                for pos in sgRNAs[gene]:
                    dist = distance_d[pos]
                    if new_zone[0] < dist < new_zone[1]:
                        VBC_score_picking = sgRNAs[gene][pos]
                        VBC_score_picking = VBC_score_picking * 0.5
                        sgRNAs[gene][pos] = VBC_score_picking
            elif hit_type == '1st':
                new_zone = [dist-0.05,dist+0.05]
                hitzone.append(new_zone)
                for pos in sgRNAs[gene]:
                    dist = distance_d[pos]
                    if new_zone[0] < dist < new_zone[1]:
                        VBC_score_picking = sgRNAs[gene][pos]
                        VBC_score_picking = VBC_score_picking * 0.95
                        sgRNAs[gene][pos] = VBC_score_picking
            # now remove pos from sgRNA[gene] and re-score all other guides
            # break after 1st guide and rerun until 6 sgRNAs picked.
            break
        #print(hitzone)