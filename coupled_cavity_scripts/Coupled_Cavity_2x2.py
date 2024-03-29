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
Ecm = 0.7294-0.00009*1j
RabiC = 0.009
XC = 0.7294-0.00009*1j

## Create empty dataframe for storing data
df_Eigen = pd.DataFrame(columns=['eigVal',
                                 'eigVec',
                                 'Coupling coeff.',
                                 'Ecm',
                                 'XC',
                                 'Ecm-XC'])

########################################
##                                    ##
##         Compute Eigenvals          ##
##                                    ##
########################################

initialC = RabiC
for a in range(1,amax+1):
    RabiC = initialC + 0.000035*(a-1)
    M = np.array([[  Ecm-0.005j,  RabiC/2    ],
                  [  RabiC/2,      XC+0.005j ]])
    eigVal, eigVec = np.linalg.eig(M)
    # append eigenvals and eigenvectors to dataframe
    df_Eigen.loc[len(df_Eigen.index)] = [eigVal, eigVec, RabiC/2, Ecm, XC, Ecm-XC]

## output the dataframe
df_Eigen.to_csv('result.csv')


########################################
##                                    ##
##          Extracting data           ##
##                                    ##
########################################
realEigVal = np.real(df_Eigen['eigVal'].values.tolist())
compEigVal = np.imag(df_Eigen['eigVal'].values.tolist())
#gainVal = np.real(df_Eigen['Gain'].values.tolist())

realEigVal1 = [vals[0] for vals in realEigVal]
realEigVal2 = [vals[1] for vals in realEigVal]

compEigVal1 = [vals[0] for vals in compEigVal]
compEigVal2 = [vals[1] for vals in compEigVal]

realCompEigVals = [[realEigVal1, realEigVal2], [compEigVal1, compEigVal2]]


########################################
##                                    ##
##            Generate plots          ##
##                                    ##
########################################

## Extract the eigenvectors
eigVecsSet = df_Eigen['eigVec']
eigVec1 = [[vecs[0][0], vecs[1][0]] for vecs in eigVecsSet]
eigVec2 = [[vecs[0][1], vecs[1][1]] for vecs in eigVecsSet]
HopfCoeffC = [vec[0]**2 for vec in eigVec1]
HopfCoeffE = [vec[1]**2 for vec in eigVec1]

Ediffs = np.real(df_Eigen['Ecm-XC'].values.tolist())
Cs = df_Eigen['Coupling coeff.']

## plot the eigenvalues
if not os.path.exists('plots'):
    os.makedirs('plots')

names = ['Eigen Val 1', 'Eigen Val 2', 'Complex', 'Real']
figs = [go.Figure(), go.Figure()]

for j in range(len(figs)):
    for i in range(len(realCompEigVals[j])):
        fig = figs[j]
        fig.add_trace(go.Scatter(x=Cs, y=realCompEigVals[j][i], mode='markers', name=names[i]))
#    if (j==0):
#        fig.add_trace(go.Scatter(x=Cs, y=df_Eigen['Ecm'], mode='lines', name='Ecm'))
#        fig.add_trace(go.Scatter(x=Cs, y=df_Eigen['XC'], mode='lines', name='XC'))
    fig.update_layout(title= names[-j-1] + ' part of Eigenvals',
                      xaxis_title='Coupling strength',
                      yaxis_title= names[-j-1] + ' part of Eigenvals',
                      showlegend=True)
    fig.write_html('plots/' + names[-j-1] + 'PartEigen.html', auto_open=True)

VecFig = go.Figure()
VecFig.add_trace(go.Scatter(x=Cs, y=HopfCoeffC, mode='markers', name='Cavity'))
VecFig.add_trace(go.Scatter(x=Cs, y=HopfCoeffE, mode='markers', name='Exciton'))
VecFig.update_layout(title= 'Hopfield Coefficients',
                  xaxis_title='Coupling strength',
                  yaxis_title='Coefficient',
                  showlegend=True)
VecFig.write_html('plots/' + '2x2_HopfieldCoeff_woGain.html', auto_open=True)
