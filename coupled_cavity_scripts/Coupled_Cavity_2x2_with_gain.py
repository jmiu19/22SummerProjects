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
Ecm1 = 1600
RabiC = 5
Ecm2 = 1600

P = np.array([[0, 1],
              [1, 0]])


## Create empty dataframe for storing data
df_Eigen = pd.DataFrame(columns=['eigVal',
                                 'eigVec',
                                 'Coupling coeff.',
                                 'Ecm1',
                                 'Ecm2',
                                 'Ecm1-Ecm2',
                                 'Gain'])

########################################
##                                    ##
##         Compute Eigenvals          ##
##                                    ##
########################################

initialG = 0
for a in range(1,amax+1):
    G = initialG + 0.01*(a-1)
    M = np.array([[  Ecm1+1j*G,    RabiC/2 ],
                  [  RabiC/2,     Ecm2-1j*G ]])
    eigVal, eigVec = np.linalg.eig(M)
    # append eigenvals and eigenvectors to dataframe
    df_Eigen.loc[len(df_Eigen.index)] = [eigVal, eigVec, RabiC/2, Ecm1, Ecm2, Ecm1-Ecm2, G]

## output the dataframe
df_Eigen.to_csv('result.csv')


########################################
##                                    ##
##          Extracting data           ##
##                                    ##
########################################
realEigVal = np.real(df_Eigen['eigVal'].values.tolist())
compEigVal = np.imag(df_Eigen['eigVal'].values.tolist())

realEigVal1 = [vals[0] for vals in realEigVal]
realEigVal2 = [vals[1] for vals in realEigVal]

compEigVal1 = [vals[0] for vals in compEigVal]
compEigVal2 = [vals[1] for vals in compEigVal]

realCompEigVals = [[realEigVal1, realEigVal2], [compEigVal1, compEigVal2]]

## Extract the eigenvectors
eigVecsSet = df_Eigen['eigVec']
eigVec1 = [[vecs[0][0], vecs[1][0]] for vecs in eigVecsSet]
eigVec2 = [[vecs[0][1], vecs[1][1]] for vecs in eigVecsSet]

HopfCoeff11 = [vec[0]*np.conj(vec[0]) for vec in eigVec1]
HopfCoeff12 = [vec[1]*np.conj(vec[1]) for vec in eigVec1]
HopfCoeff21 = [vec[0]*np.conj(vec[0]) for vec in eigVec2]
HopfCoeff22 = [vec[1]*np.conj(vec[1]) for vec in eigVec2]

orderParameter1 = []
eigenVecs = [eigVec1, eigVec2]
for i in range(len(eigVec1)):
    eta = 0
    for j in range(len(eigenVecs)):
        PTvec = np.transpose(np.matmul(P, np.conj([[eigenVecs[j][i][0]],
                                                   [eigenVecs[j][i][1]]])))[0]
        eta = eta + abs(+ np.dot(np.conj(eigenVecs[j][i]), eigenVecs[j][i])
                        - np.dot(np.conj(PTvec), eigenVecs[j][i]))
    orderParameter1.append(eta)

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
figs = [go.Figure(), go.Figure()]

for j in range(len(figs)):
    for i in range(len(realCompEigVals[j])):
        fig = figs[j]
        fig.add_trace(go.Scatter(x=Gs, y=realCompEigVals[j][i], mode='markers', name=names[i]))
    if (j==0):
        fig.add_trace(go.Scatter(x=Gs, y=df_Eigen['Ecm1'], mode='lines', name='Ecm1'))
        fig.add_trace(go.Scatter(x=Gs, y=df_Eigen['Ecm2'], mode='lines', name='Ecm2'))
    fig.update_layout(title= names[-j-1] + ' part of Eigenvals',
                      xaxis_title='Gain of cavity',
                      yaxis_title= names[-j-1] + ' part of Eigenvals',
                      showlegend=True)
    fig.write_html('plots/' + names[-j-1] + 'PartEigen.html', auto_open=True)

HopFigR = go.Figure()
HopFigR.add_trace(go.Scatter(x=Gs, y=np.real(HopfCoeff11), mode='markers', name='Mode 1 Cavity 1'))
HopFigR.add_trace(go.Scatter(x=Gs, y=np.real(HopfCoeff12), mode='markers', name='Mode 1 Cavity 2'))
HopFigR.add_trace(go.Scatter(x=Gs, y=np.real(HopfCoeff21), mode='markers', name='Mode 2 Cavity 2'))
HopFigR.add_trace(go.Scatter(x=Gs, y=np.real(HopfCoeff22), mode='markers', name='Mode 2 Cavity 2'))
HopFigR.update_layout(title= 'Hopfield Coefficients',
                  xaxis_title='Gain of cavity',
                  yaxis_title='Coefficient',
                  showlegend=True)
HopFigR.write_html('plots/' + '2x2_HopfieldCoeff.html', auto_open=True)




VecFigR = go.Figure()
VecFigR.add_trace(go.Scatter(x=Gs, y=np.real([vec[0] for vec in eigVec1]), mode='markers', name='mode 1 (1st coord.)'))
VecFigR.add_trace(go.Scatter(x=Gs, y=np.real([vec[1] for vec in eigVec1]), mode='markers', name='mode 1 (2nd coord.)'))
VecFigR.add_trace(go.Scatter(x=Gs, y=np.real([vec[0] for vec in eigVec2]), mode='markers', name='mode 2 (1st coord.)'))
VecFigR.add_trace(go.Scatter(x=Gs, y=np.real([vec[1] for vec in eigVec2]), mode='markers', name='mode 2 (2nd coord.)'))
VecFigR.update_layout(title= 'Vector Coordinates',
                  xaxis_title='Gain of cavity',
                  yaxis_title='Coefficient',
                  showlegend=True)
VecFigR.write_html('plots/' + '2x2_VectorCoordinatesReal.html', auto_open=True)

VecFigC = go.Figure()
VecFigC.add_trace(go.Scatter(x=Gs, y=np.imag([vec[0] for vec in eigVec1]), mode='markers', name='mode 1 (1st coord.)'))
VecFigC.add_trace(go.Scatter(x=Gs, y=np.imag([vec[1] for vec in eigVec1]), mode='markers', name='mode 1 (2nd coord.)'))
VecFigC.add_trace(go.Scatter(x=Gs, y=np.imag([vec[0] for vec in eigVec2]), mode='markers', name='mode 2 (1st coord.)'))
VecFigC.add_trace(go.Scatter(x=Gs, y=np.imag([vec[1] for vec in eigVec2]), mode='markers', name='mode 2 (2nd coord.)'))
VecFigC.update_layout(title= 'Vector Coordinates',
                  xaxis_title='Gain of cavity',
                  yaxis_title='Coefficient',
                  showlegend=True)
VecFigC.write_html('plots/' + '2x2_VectorCoordinatesComp.html', auto_open=True)

phaseTransPlot = go.Figure()
phaseTransPlot.add_trace(go.Scatter(x=gainVal, y=orderParameter1,
                                    mode='markers', name='1st method'))
phaseTransPlot.update_layout(title= 'Order Paramter',
                             xaxis_title='Gain of Band 1',
                             yaxis_title= 'eta',
                             showlegend=True)
phaseTransPlot.write_html('plots/' + 'orderParameter.html', auto_open=True)
