import numpy as np
import pandas as pd
import plotly.graph_objects as go
import os

## Create constant
amax = 401
C = 8 # Coupling
Ecm = 1480
Rabic = 15
RabiC = 15
XC = 1600
G=5 # Gain
L=5 # Loss

# C = input('Coupling coefficient between photonic bands: ')
# Rabic = input('Coupling coefficient between photonic band 1 and exciton: ')
# RabiC = input('Coupling coefficient between photonic band 12 and exciton: ')
# Ecm = input('Energy of each of the photonic bands: ')
# XL = input('Exciton Loss: ')


df_Eigen = pd.DataFrame(columns=['eigVal',
                                 'eigVec',
                                 'Gain',
                                 'Loss',
                                 'Coupling bt. bands',
                                 'Coupling bt. band1 and exciton',
                                 'Coupling bt. band2 and exciton',
                                 'Ecm',
                                 'XC'])

initialG = G
initialL = L
for a in range(1,amax+1):
    G = initialG + 0.1*(a-1)
    L = initialL
    M = np.array([[  Ecm+1j*G,           C,    RabiC/2],
                  [        C,     Ecm-1j*L,    Rabic/2],
                  [  RabiC/2,     Rabic/2,         XC]])
    eigVal, eigVec = np.linalg.eig(M)
    df_Eigen.loc[len(df_Eigen.index)] = [eigVal, eigVec, G, L, C, RabiC/2, Rabic/2, Ecm, XC]

df_Eigen.to_csv('result.csv')

realEigVal = np.real(df_Eigen['eigVal'].values.tolist())
compEigVal = np.imag(df_Eigen['eigVal'].values.tolist())
gainVal = np.real(df_Eigen['Gain'].values.tolist())

realEigVal1 = [vals[0] for vals in realEigVal]
realEigVal2 = [vals[1] for vals in realEigVal]
realEigVal3 = [vals[2] for vals in realEigVal]

compEigVal1 = [vals[0] for vals in compEigVal]
compEigVal2 = [vals[1] for vals in compEigVal]
compEigVal3 = [vals[2] for vals in compEigVal]




realCompEigVals = [[realEigVal1, realEigVal2, realEigVal3], [compEigVal1, compEigVal2, compEigVal3]]


RMatchV1V2Count, RMatchV1V3Count, RMatchV2V3Count = 0, 0, 0
CMatchV1V2Count, CMatchV1V3Count, CMatchV2V3Count = 0, 0, 0

RMatchV1V2Index, RMatchV1V3Index, RMatchV2V3Index = [], [], []
CMatchV1V2Index, CMatchV1V3Index, CMatchV2V3Index = [], [], []

MatchCount = [[RMatchV1V2Count, RMatchV1V3Count, RMatchV2V3Count],
              [CMatchV1V2Count, CMatchV1V3Count, CMatchV2V3Count]]
MatchIndex = [[RMatchV1V2Index, RMatchV1V3Index, RMatchV2V3Index],
              [CMatchV1V2Index, CMatchV1V3Index, CMatchV2V3Index]]

threshold = 0.5

for i in range(len(['realPart', 'complexPart'])):
    for j in range(len(realCompEigVals[i])-1):
        for t in range(len(realCompEigVals)-j):
            for k in range(len(realCompEigVals[i][j])):
                if ((realCompEigVals[i][j][k] - realCompEigVals[i][j+t+1][k]) <= threshold
                    and (realCompEigVals[i][j+t+1][k] - realCompEigVals[i][j][k]) <= threshold):
                    MatchCount[i][j*2 + 1*t] = MatchCount[i][j*2 + 1*t] + 1
                    MatchIndex[i][j*2 + 1*t].append(k)


checkIndex = ['val1 and val2', 'val1 and val3', 'val2 and val3']
excepIndex = [[], [], []]

for i in range(len(checkIndex)):
    for k in range(len(MatchIndex[0][i])):
        if any([(MatchIndex[0][i][k]<=val+1 and MatchIndex[0][i][k]>=val-1) for val in MatchIndex[1][i]]):
            excepGain = initialG + 0.1*MatchIndex[0][i][k]
            print('Found potentional exceptional point between '
                  + str(checkIndex[i])
                  + ' at around Gain = '
                  + str(excepGain))
            excepIndex[i].append(MatchIndex[0][i][k])








if not os.path.exists('plots'):
    os.makedirs('plots')

names = ['Eigen Val 1', 'Eigen Val 2', 'Eigen Val 3', 'Complex', 'Real']
figs = [go.Figure(), go.Figure()]

for j in range(len(figs)):
    for i in range(len(realCompEigVals[j])):
        fig = figs[j]
        fig.add_trace(go.Scatter(x=gainVal, y=realCompEigVals[j][i], mode='markers', name=names[i]))
        fig.add_trace(go.Scatter(x=[initialG + 0.1*(val-1) for val in excepIndex[i]],
                                 y=[realCompEigVals[j][i][t] for t in excepIndex[i]],
                                 mode='markers', name='Exceptional Point', marker=dict(size=19, color='black')))
        fig.update_layout(title= names[-j-1] + ' part of Eigenvals',
                          xaxis_title='Gain of Band 1',
                          yaxis_title= names[-j-1] + ' part of Eigenvals',
                          showlegend=True)
    fig.write_html('plots/' + names[-j-1] + 'PartEigen.html', auto_open=True)
