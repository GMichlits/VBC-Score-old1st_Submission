__author__ = 'georg.michlits'

import numpy as np
import matplotlib.pyplot as plt


name = 'V2_hap.insuf._red hap.suf_blue 60-100%KO'
hap_KO_prob = 0.80

model_noise_out = open('model_60-100LOF.txt','w')

Noise_lv = 0.3

hap_KO_prob = 0.6+(np.random.random())*0.4
dip_suff = hap_KO_prob**2
dip_insuff = hap_KO_prob+(1-hap_KO_prob)*hap_KO_prob

model_noise_out.write('X_recessive' + '\t' + 'Y_recessive' + '\t' + 'X_dominant' + '\t' + 'Y_dominant' + '\t' + 'X_100LOF' + '\t' + 'Y_100LOF')
x_ls_suf = []
y_ls_suf = []
x_ls_insuf = []
y_ls_insuf = []
x_ls_bl = []
y_ls_bl = []
for i in range(20000):
    ff_rand = ((np.random.random()+np.random.random()+np.random.random()+np.random.random()+np.random.random())/5)
    if np.random.randint(10) == 5:
        ff_rand = ((np.random.random())/np.random.randint(1,3))
    if np.random.randint(600) == 52:
        ff_rand = ((np.random.random())/np.random.randint(1,50))
    if np.random.randint(1000) == 52:
        ff_rand = ((np.random.random())/np.random.randint(1,200))
    if np.random.randint(2000) == 52:
        ff_rand = ((np.random.random())/np.random.randint(1,500))
    if np.random.randint(3500) == 52:
        ff_rand = ((np.random.random())/np.random.randint(1,2000))
    LFC = np.log2(ff_rand/0.5)
    hap_KO_prob =0.6+np.random.random()*0.4
    dip_suff = hap_KO_prob**2
    dip_insuff = hap_KO_prob+(1-hap_KO_prob)*hap_KO_prob
    y_ns = (np.random.rand()+np.random.rand()+np.random.rand())/3*Noise_lv
    x_ns = (np.random.rand()+np.random.rand()+np.random.rand())/3*Noise_lv
    x_suf = np.log2((0.5-(0.5-ff_rand)*dip_suff)/0.5)+y_ns
    y_suf = np.log2((0.5-(0.5-ff_rand)*hap_KO_prob)/0.5)+x_ns
    x_insuf = np.log2((0.5-(0.5-ff_rand)*dip_insuff)/0.5)+y_ns
    y_insuf = y_suf
    x_bl = LFC+x_ns
    y_bl = LFC+y_ns
    y_ls_suf.append(y_suf)
    x_ls_suf.append(x_suf)
    y_ls_insuf.append(y_insuf)
    x_ls_insuf.append(x_insuf)
    x_ls_bl.append(x_bl)
    y_ls_bl.append(y_bl)
    model_noise_out.write('\n' + str(x_suf) + '\t' + str(y_suf) + '\t' + str(x_insuf) + '\t' + str(y_insuf) + '\t' + str(x_bl) + '\t' + str(y_bl))

plt.figure(figsize=(12,10))
plt.scatter(x_ls_bl,y_ls_bl, c='grey', alpha=0.2, s=120)
plt.scatter(x_ls_suf,y_ls_suf, c='green', alpha=0.2, s=120)
plt.scatter(x_ls_insuf,y_ls_suf, c='purple', alpha=0.2, s=120)
plt.xlim(-8,2)
plt.tick_params(labelsize=30)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.ylim(-8,2)
#plt.show()
plt.savefig(name + 'purplegreen_v2.png')
plt.savefig(name + 'purplegreen_v2.pdf')