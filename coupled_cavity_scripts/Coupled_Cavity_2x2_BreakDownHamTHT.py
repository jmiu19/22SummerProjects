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

amax = 800 # number of iteration when computing eigenvalues

## Default constants
C = 8 # Coupling
detuning = 7.5
XC = 1650
RabiC = 10
G= 0 # Gain
Gu = 30
L='G' # Loss

# # ask user for constants input
# C = float(input('Coupling coefficient between photonic cavities: '))
# RabiC = float(input('Coupling coefficient between photonic cavities and exciton: '))
# detuning = float(input('Detuning of the photonic cavities: '))
# XC = float(input('Energy of the exciton: '))
# G = float(input('Lowest gain: '))
# Gu = float(input('Highest gain: '))
# L = input('Loss of the cavity: (input G if setting loss = gain)   ')

Ecm = XC + detuning

## Create empty dataframe for storing data
df_Eigen = pd.DataFrame(columns=['eigVal1', 'eigVec1',
                                 'eigVal2', 'eigVec2',
                                 'eigVal3', 'eigVec3',
                                 'RabiC/2', 'Xc',
                                 'Ecm', 'Gain', 'Coupling'])


########################################
##                                    ##
##         Compute Eigenvals          ##
##                                    ##
########################################

initialG = 0
for a in range(1,amax+1):
    G = initialG + (Gu-initialG)/(amax)*(a)
    M1 = np.array([[ Ecm+C,     1j*G ],
                   [  1j*G,    Ecm-C ]])
    M2 = np.array([[              Ecm+C,  RabiC/(np.sqrt(2)) ],
                   [ RabiC/(np.sqrt(2)),                  XC ]])
    M3 = np.array([[ Ecm-C,     0 ],
                   [     0,    XC ]])
    eigVal1, eigVec1 = np.linalg.eig(M1)
    eigVal2, eigVec2 = np.linalg.eig(M2)
    eigVal3, eigVec3 = np.linalg.eig(M3)
    # append eigenvals and eigenvectors to dataframe
    df_Eigen.loc[len(df_Eigen.index)] = [eigVal1, eigVec1,
                                         eigVal2, eigVec2,
                                         eigVal3, eigVec3,
                                         RabiC/2,      XC,
                                         Ecm, G, C]

## output the dataframe
df_Eigen.to_csv('result.csv')


########################################
##                                    ##
##          Extracting data           ##
##                                    ##
########################################

# a list that contains eigenvals for all 3 systems
EigenVals = []

for i in range(len(['M1', 'M2', 'M3'])):
    realEigVal = np.real(df_Eigen['eigVal'+str(i+1)].values.tolist())
    compEigVal = np.imag(df_Eigen['eigVal'+str(i+1)].values.tolist())

    realEigValComp1 = [vals[0] for vals in realEigVal]
    realEigValComp2 = [vals[1] for vals in realEigVal]
    compEigValComp1 = [vals[0] for vals in compEigVal]
    compEigValComp2 = [vals[1] for vals in compEigVal]

    realCompEigVals = [[realEigValComp1, realEigValComp2],
                       [compEigValComp1, compEigValComp2]]
    EigenVals.append(realCompEigVals)


## Extract the eigenvectors
# a list that contains eigenvals for all 3 systems
EigenVecs = []
HopfCoeff = []

for i in range(len(['M1', 'M2', 'M3'])):
    eigVecsSet = df_Eigen['eigVec'+str(i+1)]
    eigVec1 = [[vecs[0][0], vecs[1][0]] for vecs in eigVecsSet]
    eigVec2 = [[vecs[0][1], vecs[1][1]] for vecs in eigVecsSet]
    EigenVecs.append([eigVec1, eigVec2])

    HopfCoeff11 = [vec[0]*np.conj(vec[0]) for vec in eigVec1]
    HopfCoeff12 = [vec[1]*np.conj(vec[1]) for vec in eigVec1]
    HopfCoeff21 = [vec[0]*np.conj(vec[0]) for vec in eigVec2]
    HopfCoeff22 = [vec[1]*np.conj(vec[1]) for vec in eigVec2]
    HopfCoeff.append([HopfCoeff11, HopfCoeff12, HopfCoeff21, HopfCoeff22])



Gs = np.real(df_Eigen['Gain'].values.tolist())


########################################
##                                    ##
##            Generate plots          ##
##                                    ##
########################################

## plot the eigenvalues
if not os.path.exists('plots'):
    os.makedirs('plots')

names = ['Eigen Val 1', 'Eigen Val 2', 'Complex', 'Real']
figs = [[go.Figure(), go.Figure()],
        [go.Figure(), go.Figure()],
        [go.Figure(), go.Figure()]]

for n in range(len(EigenVals)):
    realCompEigVals = EigenVals[n]
    for j in range(len(EigenVals[n])):
        for i in range(len(realCompEigVals[j])):
            fig = figs[n][j]
            fig.add_trace(go.Scatter(x=Gs, y=realCompEigVals[j][i], mode='markers', name=names[i]))
        if (j==0):
            fig.add_trace(go.Scatter(x=Gs, y=df_Eigen['Ecm'], mode='lines', name='Ecm'))
        fig.update_layout(title= names[-j-1] + ' part of Eigenvals for system M'+str(n+1),
                          xaxis_title='Gain of cavity',
                          yaxis_title= names[-j-1] + ' part of Eigenvals',
                          showlegend=True)
        fig.write_html('plots/' + names[-j-1] + 'PartEigenM' + str(n+1) + '.html', auto_open=True)

HopFig = [go.Figure(), go.Figure(), go.Figure()]
HopFigNames = ['Mode 1 Cavity 1', 'Mode 1 Cavity 2', 'Mode 2 Cavity 2', 'Mode 2 Cavity 2']
for n in range(len(['M1', 'M2', 'M3'])):
    fig = HopFig[n]
    for i in range(len(HopfCoeff[n])):
        fig.add_trace(go.Scatter(x=Gs, y=np.real(HopfCoeff[n][i]),
                                 mode='markers', name=HopFigNames[i]))
    fig.update_layout(title= 'Hopfield Coefficients for system M' + str(n+1),
                      xaxis_title='Gain of cavity',
                      yaxis_title='Coefficient',
                      showlegend=True)
    fig.write_html('plots/' + '2x2_HopfieldCoeffM' + str(n+1) + '.html', auto_open=True)
