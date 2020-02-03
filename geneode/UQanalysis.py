###------Network Simulator------###
###------Shulin Cao------###
###------CMRG, UC San Diego------###

###import packages###
import pandas as pd
import collections
import timeit
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import networkx as nx
import matplotlib.mlab as mlab
import statistics
import random
import numpy.linalg
import numpy as np
import sys
from scipy.optimize import minimize
elapsed_time = timeit.default_timer()
from sklearn.linear_model import LinearRegression
from sklearn import cluster
import seaborn as sns
sns.set()
from sklearn import datasets
from sklearn.metrics import r2_score
from matplotlib import pylab
from matplotlib import font_manager
import chaospy as cp 
import uncertainpy as un
def get_reactors(reac):
    reac_split = reac.split(' ')
    reactors = []
    for k in reac_split:
        if k != '&' and k!= '=>':
            reactors.append(k)
    return reactors[:-1]


def Hill(reactor, n, EC50):
    B = (EC50**n-1)/(2*EC50**n-1)
    C = (B-1)**(1/n)
    if reactor[0] == '!':
        return (1-B*globals()['{}'.format(reactor[1:])]**n/(C**n + globals()['{}'.format(reactor[1:])]**n))
    else:
        return B*globals()['{}'.format(reactor)]**n/(C**n + globals()['{}'.format(reactor)]**n)


def OR(reaction_list):
    tera = (-1)**(len(reaction_list)+1)
    for k in reaction_list:
        weight, n, EC50 = reaction_list[k]
        final = weight
        for j in get_reactors(k):
            final *= Hill(j, n, EC50)
        tera *= (final-1)
    tera +=1
    return tera
    
def inte(state, t, reaction_dict):
    for i in range(len(node_ID)):
        globals()['{}'.format(node_ID[i])] = state[i]
    for i in range(len(node_ID)):
        if len(reaction_dict[node_ID[i]]) == 1:
            reactors = get_reactors(list(reaction_dict[node_ID[i]].keys())[0])
            weight, n, EC50 = reaction_dict[node_ID[i]][list(reaction_dict[node_ID[i]].keys())[0]]
            TF = 1
            for j in reactors:
                TF *= Hill(j, n, EC50)
            globals()['{}'.format(node_ID[i] + 'd')] = (TF*weight*Ymax[i]-globals()['{}'.format(node_ID[i])])/tau[i]
        else:
            TF = OR(reaction_dict[node_ID[i]])
            globals()['{}'.format(node_ID[i] + 'd')] = (TF*Ymax[i]-globals()['{}'.format(node_ID[i])])/tau[i]
    return [globals()['{}'.format(k+ 'd')] for k in node_ID]
    
reactions_raw = pd.read_excel('/home/shulincao/Desktop/projects/demo/MTv29-philip-20170921-for-network.xlsx', sheet_name = 1, skiprows = 1, header = 0)
species_raw = pd.read_excel('/home/shulincao/Desktop/projects/demo/MTv29-philip-20170921-for-network.xlsx', sheet_name = 0, skiprows = 1, header = 0)
pmid = reactions_raw['PMID'].tolist()

species = species_raw[['ID', 'Yinit', 'Ymax', 'tau']]
node_ID = species['ID'].tolist()
Yinit = species['Yinit'].tolist()
Ymax = species['Ymax'].tolist()
tau = species['tau'].tolist()

species_dict = dict()
for k in range(len(species)):
    #lis = species.loc[k, ['Yinit', 'Ymax', 'tau']].tolist()
    species_dict[species.loc[k, 'ID']] = species.loc[k, ['Yinit', 'Ymax', 'tau']].tolist()

state0 = []
for k in range(len(node_ID)):
    state0.append(Yinit[k])  #solve_ivp
    
t = np.arange(0.0, 60.0*4, 0.1)

validation = {}
"""
validation['Stretch'] = {
    'aMHC':-1,
    'ANP':1,
    'Ao':1,
    'bMHC':1,
    'BNP':1,
    'CellArea':1,
    'PrSynth':1,
    'sACT':1,
    'SERCA':-1,
    'Akt':1,
    'AngII':1,
    'AP1':1,
    'Ca':1,
    'CaN':1,
    'cFos':1,
    'cJun':1,
    'cMyc':1,
    'CREB':1,
    'Cx43':1,
    'DAG':1,
    'EGFR':1,
    'ERK12':1,
    'FAK':1,
    'FHL1':1,
    'GATA4':1,
    'gp130':1,
    'GSK3b':-1,
    'IP3':1,
    'JAK':1,
    'JNK':1,
    'Lmcd1':1,
    'MEF2':1,
    'MEK12':1,
    'MLP':1,
    'MRTF':1,
    'mTor':1,
    'MuRF':1,
    'NFAT':1,
    'NFkB':1,
    'NOS':1,
    'p38':1,
    'p70s6k':1,
    'PI3K':1,
    'PKC':1,
    'Rac1':1,
    'Raf1':1,
    'Ras':1,
    'RhoA':1,
    'RhoGEF':1,
    'Src':1,
    'SRF':1,
    'STAT':1
}
"""
validation['Akt'] = {'ERK12':-1}
validation['AP1'] = {'BNP':0}
validation['AT1R'] = {
    'ANP':-1,
    'Ao':-1,
    'BNP':-1,
    'CellArea':-1,
    'cFos':-1,
    'cJun':0,
    'Cx43':-1,
    'ERK12':-1,
    'JNK':1,
    'Raf1':-1,
    'sACT':-1,
    'STAT':-1
}
validation['Ca'] = {
    'cFos':-1,
    'cJun':0,
    'STAT':-1
}
validation['CaN'] = {
    'ANP':-1
}
validation['EGFR'] = {
    'BNP':-1,
    'ERK12':-1,
    'JNK':0,
    'MEK12':-1,
    'Ras':-1
}
validation['ET1R'] = {
    'ANP':-1,
    'BNP':-1,
    'cFos':-1,
    'STAT':0
}
validation['FAK'] = {
    'Akt':-1,
    'ANP':-1,
    'bMHC':-1,
    'CellArea':-1,
    'cJun':-1,
    'cMyc':-1,
    'ERK12':-1,
    'JNK':1,
    'MEF2':-1,
    'mTor':-1,
    'p70s6k':-1,
    'Src':-1,
}
validation['Ga1213'] = {
    'RhoA':-1,
    'RhoGEF':-1
}
validation['GATA4'] = {
    'BNP':-1
}
validation['gp130'] = {
    'STAT':-1
}
validation['Integrin'] = {
    'ERK12':-1,
    'FAK':-1,
    'JNK':-1,
    'p38':-1,
    'RhoA':-1,
    'RhoGEF':-1
}
validation['JAK'] = {
    'STAT':-1
}
validation['JNK'] = {
    'ANP':-1,
    'Ao':1,
    'cJun':-1,
    'ERK12':-1
}
validation['Lmcd1'] = {
    'CellArea':-1
}
validation['LTCC'] = {
    'aMHC':-1,
    'ANP':-1,
    'bMHC':-1,
    'Ca':-1,
    'CaN':-1,
    'PrSynth':-1,
    'SERCA':0
}
validation['MEK12'] = {
    'BNP':-1,
    'Cx43':-1,
    'ERK12':-1
}
validation['MLP'] = {
    'BNP':-1,
    'NFAT':-1,
    'PrSynth':-1
}
validation['MRTF'] = {
    'bMHC':-1,
    'BNP':-1
}
validation['NCX'] = {
    'ANP':-1,
    'CaN':-1,
    'PrSynth':-1
}
validation['NHE'] = {
    'ANP':-1,
    'CaN':-1,
    'ERK12':-1,
    'PrSynth':-1,
    'Raf1':-1,
    'STAT':-1
}
validation['p38'] = {
    'Ao':-1,
    'PrSynth':-1
}
validation['PI3K'] = {
    'Akt':-1,
    'BNP':-1,
    'ERK12':-1,
    'JNK':0,
    'NOS':-1,
    'Ras':-1
}
validation['PKC'] = {
    'cFos':-1,
    'Cx43':0,
    'ERK12':-1,
    'Raf1':-1,
    'STAT':-1
}
validation['PLC'] = {
    'Ca':-1,
    'cFos':-1,
    'IP3':-1
}
validation['Rac1'] = {
    'ERK12':-1
}
validation['Raf1'] = {
    'ERK12':-1
}
validation['Ras'] = {
    'ERK12':0,
    'JNK':0,
    'MEK12':-1,
    'p38':-1
}
validation['RhoGEF'] = {
    'ANP':-1,
    'bMHC':-1,
    'CellArea':-1,
    'MRTF':-1,
    'RhoA':-1
}
validation['RhoA'] = {
    'Akt':-1,
    'ANP':-1,
    'bMHC':-1,
    'BNP':-1,
    'cFos':-1,
    'ERK12':-1,
    'FAK':-1,
    'MRTF':-1,
    'PrSynth':-1,
    'sACT':-1
}
validation['Src'] = {
    'ANP':-1,
    'FAK':-1,
    'p38':-1
}
validation['Titin'] = {
    'MuRF':1
}

def acc_test(simu_data, reaction_dict, threshold):
    s = 0
    s_all = 0
    s_stretch = 0
    s_all_stretch = 0
    for i in validation:
        if i == 'Stretch':
            for j in validation[i]:
                s_all += 1
                if i == 'Stretch':
                    s_all_stretch += 1
                if simu_data[2399, node_ID.index(j)] - simu_data[0, node_ID.index(j)] >= threshold \
                    and validation[i][j]==1:
                    s += 1
                    s_stretch += 1
                elif simu_data[2399, node_ID.index(j)] - simu_data[0, node_ID.index(j)] <= -threshold \
                    and validation[i][j]==-1:
                    s += 1
                    s_stretch += 1
                elif simu_data[2399, node_ID.index(j)] - simu_data[0, node_ID.index(j)] >= -threshold \
                    and simu_data[2399, node_ID.index(j)] - simu_data[0, node_ID.index(j)] <= threshold \
                    and validation[i][j]==0:
                    s += 1
                    s_stretch += 1
        else:
            Ymax[node_ID.index(i)] = 0
            simu_data_change = odeint(inte, state0, t, args = (reaction_dict, ))
            Ymax[node_ID.index(i)] = 1
            for j in validation[i]:
                s_all += 1
                if simu_data_change[2399, node_ID.index(j)] - simu_data[2399, node_ID.index(j)] >= threshold \
                    and validation[i][j]==1:
                    s += 1
                elif simu_data_change[2399, node_ID.index(j)] - simu_data[2399, node_ID.index(j)] <= -threshold \
                    and validation[i][j]==-1:
                    s += 1
                elif simu_data_change[2399, node_ID.index(j)] - simu_data[2399, node_ID.index(j)] >= -threshold \
                    and simu_data_change[2399, node_ID.index(j)] - simu_data[2399, node_ID.index(j)] <= threshold \
                    and validation[i][j]==0:
                    s += 1
    return s/s_all
    
    def acc_stretch_val(val0, val1, val2, val3, val4, val5, val6, val7, val8, val9):
    vali_series = [val0, val1, val2, val3, val4, val5, val6, val7, val8, val9]
    reactions = {
        'rule':reactions_raw['rule'].tolist(), 
        'weight':reactions_raw['weight'].tolist(), 
        'n':reactions_raw['n'].tolist(), 
        'EC50':reactions_raw['EC50'].tolist(), 
    }
    val0 = 0 if val0<=0.6 else 1
    val1 = 0 if val1<=0.6 else 1
    val2 = 0 if val2<=0.6 else 1
    val3 = 0 if val3<=0.6 else 1
    val4 = 0 if val4<=0.6 else 1
    val5 = 0 if val5<=0.6 else 1
    val6 = 0 if val6<=0.6 else 1
    val7 = 0 if val7<=0.6 else 1
    val8 = 0 if val8<=0.6 else 1
    val0 = 0 if val0<=0.6 else 1
    validation['AP1']['BNP'] = val0
    validation['AT1R']['cJun'] = val1
    validation['Ca']['cJun'] = val2
    validation['EGFR']['JNK'] = val3
    validation['ET1R']['STAT'] = val4
    validation['LTCC']['SERCA'] = val5
    validation['PI3K']['JNK'] = val6
    validation['PKC']['Cx43'] = val7
    validation['Ras']['ERK12'] = val8
    validation['Ras']['JNK'] = val9   
    reactions = pd.DataFrame(data = reactions)
    reaction_dict = collections.defaultdict(dict)
    for k in range(len(reactions)):
        node = reactions.loc[k, 'rule'].split(' ')
        reaction_dict[node[-1]][reactions.loc[k, 'rule']] = reactions.loc[k, ['weight', 'n', 'EC50']].tolist() 
    reaction_dict['Stretch']['=> Stretch'] = [0.7, 1.4, 0.5]
    simu_data = odeint(inte, state0, t, args = (reaction_dict, ))
    accuracy = acc_test(simu_data, reaction_dict, threshold)
    vali_series.append(accuracy)
    info = {'parameters_func':vali_series}
    return t, accuracy, info

def parameters_func(t, accuracy, info):
    return t, info['parameters_func']

feature_list = [parameters_func]

threshold = 0.05

accuracy_sample_mc_nc = {}
parameters_nc = {"val"+str(i): cp.Uniform(0, 1) for i in range(10)}
uqmodel = un.Model(run=acc_stretch_val, labels=["Time (min)", 'Accuracy'])
UQ = un.UncertaintyQuantification(
    model=uqmodel, 
    parameters=parameters_nc, 
    features=feature_list
)

mc_samples = [100,200]
for sample in mc_samples:
    data_sample_mc = UQ.quantify(method = 'mc', nr_mc_samples = sample, plot = None)
    current_res = []
    a = data_sample_mc['parameters_func']['evaluations']
    for i in range(len(a)):
        if a[i][-1] > 0:
            current_res.append(a[i])
    accuracy_sample_mc_nc[sample] = current_res
    
val = [[] for i in range(10)]
acc_val = []
val_sample = []
for i in accuracy_sample_mc_nc:
    for j in accuracy_sample_mc_nc[i]:
        for m in range(10):
            val[m].append(j[m])
        val_sample.append(len(accuracy_sample_mc_nc[i]))
        acc_val.append(j[-1])
val_result = {'val'+str(i): val[i] for i in range(10)}
val_result['Sample'] = val_sample
val_result['Accuracy'] = acc_val
val_result = pd.DataFrame(val_result)

val_result.to_csv('val_nc_power06_montecarlo.csv', index = False)
