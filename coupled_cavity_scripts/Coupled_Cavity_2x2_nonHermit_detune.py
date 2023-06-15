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
E2 = 0
detuning = -10
E1 = E2 + detuning
C = 1
r=5

## Create empty dataframe for storing data
df_Eigen = pd.DataFrame(columns=['eigVal',
                                 'eigVec',
                                 'Coupling coeff.',
                                 'E1',
                                 'E2',
                                 'E1-E2'])

########################################
##                                    ##
##         Compute Eigenvals          ##
##                                    ##
########################################

initialE1 = E1
for a in range(1,amax+1):
    E1 = initialE1 + 0.1*(a-1)
    M = np.array([[   E1,   C ],
                  [  r*C,   E2]])
    eigVal, eigVec = np.linalg.eig(M)
    # append eigenvals and eigenvectors to dataframe
    df_Eigen.loc[len(df_Eigen.index)] = [eigVal, eigVec, C, E1, E2, E1-E2]

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

Ediffs = np.real(df_Eigen['E1-E2'].values.tolist())


## plot the eigenvalues
if not os.path.exists('plots'):
    os.makedirs('plots')

names = ['Eigen Val 1', 'Eigen Val 2', 'Complex', 'Real']
figs = [go.Figure(), go.Figure()]

for j in range(len(figs)):
    for i in range(len(realCompEigVals[j])):
        fig = figs[j]
        fig.add_trace(go.Scatter(x=Ediffs, y=realCompEigVals[j][i], mode='markers', name=names[i]))
    if (j==0):
        fig.add_trace(go.Scatter(x=Ediffs, y=df_Eigen['E1'], mode='lines', name='E1'))
        fig.add_trace(go.Scatter(x=Ediffs, y=df_Eigen['E2'], mode='lines', name='E2'))
    fig.update_layout(title= names[-j-1] + ' part of Eigenvals',
                      xaxis_title='Energy difference between eE2iton and cavity',
                      yaxis_title= names[-j-1] + ' part of Eigenvals',
                      showlegend=True)
    fig.write_html('plots/' + names[-j-1] + 'PartEigen.html', auto_open=True)

VecFig = go.Figure()
VecFig.add_trace(go.Scatter(x=Ediffs, y=HopfCoeffC, mode='markers', name='Cavity'))
VecFig.add_trace(go.Scatter(x=Ediffs, y=HopfCoeffE, mode='markers', name='EE2iton'))
VecFig.update_layout(title= 'Hopfield Coefficients',
                  xaxis_title='Energy difference between eE2iton and cavity',
                  yaxis_title='Coefficient',
                  showlegend=True)
VecFig.write_html('plots/' + '2x2_HopfieldCoeff_woGain.html', auto_open=True)
