import numpy as np
import pandas as pd
import plotly.graph_objects as go
import os

########################################
##                                    ##
##       Setting up the model         ##
##                                    ##
########################################

## Create constant

amax = 401 # number of iteration when computing eigenvalues

## Default constants
C = 8 # Coupling
Ecm = 1480
Rabic = 15
RabiC = 15
XC = 1600
G=5 # Gain
L=5 # Loss

## ask user for constants input
# C = input('Coupling coefficient between photonic cavities: ')
# Rabic = input('Coupling coefficient between photonic cavity 1 and exciton: ')
# RabiC = input('Coupling coefficient between photonic cavity 2 and exciton: ')
# Ecm = input('Energy of each of the photonic cavities: ')
# XC = input('Energy of the exciton: ')
# G = input('Lowest gain: ')
# L = input('Loss of the cavity: ')

## Create empty dataframe for storing data
df_Eigen = pd.DataFrame(columns=['eigVal',
                                 'eigVec',
                                 'Gain',
                                 'Loss',
                                 'Coupling bt. bands',
                                 'Coupling bt. band1 and exciton',
                                 'Coupling bt. band2 and exciton',
                                 'Ecm',
                                 'XC'])

########################################
##                                    ##
##         Compute Eigenvals          ##
##                                    ##
########################################

initialG = G
initialL = L
for a in range(1,amax+1):
    G = initialG + 0.1*(a-1)
    L = initialL
    M = np.array([[  Ecm+1j*G,           C,    RabiC/2],
                  [        C,     Ecm-1j*L,    Rabic/2],
                  [  RabiC/2,     Rabic/2,         XC]])
    eigVal, eigVec = np.linalg.eig(M)
    # append eigenvals and eigenvectors to dataframe
    df_Eigen.loc[len(df_Eigen.index)] = [eigVal, eigVec, G, L, C, RabiC/2, Rabic/2, Ecm, XC]

## output the dataframe
df_Eigen.to_csv('result.csv')


########################################
##                                    ##
##          Extracting data           ##
##                                    ##
########################################
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


########################################
##                                    ##
##       Find exceptional points      ##
##                                    ##
########################################

## since the results are computed numerically
## two identical eigenvalues might have a very small numerical difference
## we define the maximal difference two eigenvals can have to be considered as being equal
threshold = 0.5

## we check the real part and the complex part separately
# for each (real and complex) part, we iterate all possible pairs of real/complex parts of eigenvals
# to check if any two real/complex parts of eigenvals are the same
for i in range(len(['realPart', 'complexPart'])):
    for j in range(len(realCompEigVals[i])-1):
        for t in range(len(realCompEigVals)-j):
            for k in range(len(realCompEigVals[i][j])):
                ## check if two real/complex parts of eigenvals are the same, if true, output the index of that part
                ## need (-threshold <= val1 - val2 <= threshold)
                if ((realCompEigVals[i][j][k] - realCompEigVals[i][j+t+1][k]) <= threshold
                    and (realCompEigVals[i][j+t+1][k] - realCompEigVals[i][j][k]) <= threshold):
                    MatchCount[i][j*2 + 1*t] = MatchCount[i][j*2 + 1*t] + 1
                    MatchIndex[i][j*2 + 1*t].append(k)


checkIndex = ['val1 and val2', 'val1 and val3', 'val2 and val3']
excepIndex = [[], [], []]

## Check the (matched complex/real) part indices to see if there are pairs of eigenval that have identical real and complex parts
for i in range(len(checkIndex)):
    for k in range(len(MatchIndex[0][i])):
        if any([(MatchIndex[0][i][k]<=val+1 and MatchIndex[0][i][k]>=val-1) for val in MatchIndex[1][i]]):
            excepGain = initialG + 0.1*MatchIndex[0][i][k]
            print('Found potentional exceptional point between '
                  + str(checkIndex[i])
                  + ' at around Gain = '
                  + str(excepGain))
            excepIndex[i].append(MatchIndex[0][i][k])


########################################
##                                    ##
##            Generate plots          ##
##                                    ##
########################################

## Extract the eigenvectors
eigVecsSet = df_Eigen['eigVec']
gainVal = np.real(df_Eigen['Gain'].values.tolist())
eigVec1, eigVec2, eigVec3 = [vecs[0] for vecs in eigVecsSet], [vecs[1] for vecs in eigVecsSet], [vecs[2] for vecs in eigVecsSet]


################# 3D plots ######################
# ## plot the real part of the eigenvectors     #
# vecFig = go.Figure()                          #
# vecFig.add_trace(go.Cone(                     #
#     x=[np.real(vec[0]) for vec in eigVec1],   #
#     y=[np.real(vec[1]) for vec in eigVec1],   #
#     z=[np.real(vec[2]) for vec in eigVec1],   #
#     u=[np.real(vec[0]) for vec in eigVec1],   #
#     v=[np.real(vec[1]) for vec in eigVec1],   #
#     w=[np.real(vec[2]) for vec in eigVec1],   #
#     sizeref=1000,                             #
#     anchor="tip",                             #
#     name="Eigvec 1",                          #
#     sizemode='absolute'))                     #
#                                               #
# vecFig.add_trace(go.Cone(                     #
#     x=[np.real(vec[0]) for vec in eigVec2],   #
#     y=[np.real(vec[1]) for vec in eigVec2],   #
#     z=[np.real(vec[2]) for vec in eigVec2],   #
#     u=[np.real(vec[0]) for vec in eigVec2],   #
#     v=[np.real(vec[1]) for vec in eigVec2],   #
#     w=[np.real(vec[2]) for vec in eigVec2],   #
#     sizeref=1000,                             #
#     anchor="tip",                             #
#     name="Eigvec 2",                          #
#     sizemode='absolute'))                     #
#                                               #
# vecFig.add_trace(go.Cone(                     #
#     x=[np.real(vec[0]) for vec in eigVec3],   #
#     y=[np.real(vec[1]) for vec in eigVec3],   #
#     z=[np.real(vec[2]) for vec in eigVec3],   #
#     u=[np.real(vec[0]) for vec in eigVec3],   #
#     v=[np.real(vec[1]) for vec in eigVec3],   #
#     w=[np.real(vec[2]) for vec in eigVec3],   #
#     sizeref=5000,                             #
#     anchor="tip",                             #
#     name="Eigvec 3",                          #
#     sizemode='absolute'))                     #
#                                               #
#vecFig.update_traces(hoverinfo="u+v+w+name",   #
#                     showscale=False,)         #
#                                               #
#vecFig.write_html('plots/Eigenvec.html', auto_open=True)
#################################################


## plot the eigenvalues
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
