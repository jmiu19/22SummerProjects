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
C = 10 # Coupling
detuning = -50
h_detuning = 40000
XC = 0
Rabic = 10
RabiC = 10
r = 5

# # ask user for constants input
# C = float(input('value of C: '))
# Rabic = float(input('value of 2*Omega: '))
# RabiC = Rabic
# XC = float(input('energy of the exciton: '))
# detuning = float(input('lowest detuning: '))
# h_detuning = float(input('highest detuning: '))
# r = float(input('ratio of couplings'))

Ecm = XC + detuning

## Create empty dataframe for storing data
df_Eigen = pd.DataFrame(columns=['eigVal',
                                 'eigVec',
                                 'Coupling bt. bands',
                                 'Coupling bt. bands and exciton',
                                 'Ecm',
                                 'XC',
                                 'detuning'])


########################################
##                                    ##
##          Define Operator           ##
##                                    ##
########################################

P = np.array([[0, 1, 0],
              [1, 0, 0],
              [0, 0, 1]])



########################################
##                                    ##
##         Compute Eigenvals          ##
##                                    ##
########################################

initialDel = detuning

# Compute eigenvals for H while sweeping detuning
for a in range(0,amax+1):
    detuning = initialDel + (h_detuning-initialDel)/(amax)*(a)
    Ecm = XC + detuning
    H = np.array([[       Ecm,          r*C,    RabiC/2],
                  [         C,          Ecm,    Rabic/2],
                  [   RabiC/2,      Rabic/2,         XC]])
    eigVal, eigVec = np.linalg.eig(H)
    # append eigenvals and eigenvectors to dataframe
    df_Eigen.loc[len(df_Eigen.index)] = [eigVal, eigVec, C, Rabic/2, Ecm, XC, detuning]

##### NOTE #####
# the output eigenvectors is a 3x3 numpy array matrix, each column is an eigenvector,
# but the numpy matrix itself is an array, or a list of rows
# that is eigVec[0] gives the first row, not the first column

## output the dataframe
df_Eigen.to_csv('result.csv')


########################################
##                                    ##
##          Extracting data           ##
##                                    ##
########################################
realEigVal = np.real(df_Eigen['eigVal'].values.tolist())
compEigVal = np.imag(df_Eigen['eigVal'].values.tolist())
DelVal = np.real(df_Eigen['detuning'].values.tolist())
## Extract the eigenvectors
eigVecsSet = df_Eigen['eigVec']
DelVal = np.real(df_Eigen['detuning'].values.tolist())

##### NOTE ##### refers to notes in the Compute Eigenvals section
# the output eigenvectors is a 3x3 numpy array matrix, each column is an eigenvector,
# but the numpy matrix itself is an array, or a list of rows
# that is eigVec[0] gives the first row, not the first column
eigVec1 = [[vecs[0][0], vecs[1][0], vecs[2][0]] for vecs in eigVecsSet]
eigVec2 = [[vecs[0][1], vecs[1][1], vecs[2][1]] for vecs in eigVecsSet]
eigVec3 = [[vecs[0][2], vecs[1][2], vecs[2][2]] for vecs in eigVecsSet]

## sort eigenvals and eigenvecs based on real part of eigenvals
for i in range(len(realEigVal)):
    for list in [realEigVal, compEigVal, eigVec1, eigVec2, eigVec3]:
        list[i] = [x for _, x in sorted(zip(realEigVal[i], list[i]))]

realEigVal1 = [vals[0] for vals in realEigVal]
realEigVal2 = [vals[1] for vals in realEigVal]
realEigVal3 = [vals[2] for vals in realEigVal]
compEigVal1 = [vals[0] for vals in compEigVal]
compEigVal2 = [vals[1] for vals in compEigVal]
compEigVal3 = [vals[2] for vals in compEigVal]
realCompEigVals = [[realEigVal1, realEigVal2, realEigVal3], [compEigVal1, compEigVal2, compEigVal3]]

## Compute the hopfield coefficients
hopfieldCoeff11 = np.real([component*np.conj(component) for component in [vec[0] for vec in eigVec1]])
hopfieldCoeff12 = np.real([component*np.conj(component) for component in [vec[1] for vec in eigVec1]])
hopfieldCoeff13 = np.real([component*np.conj(component) for component in [vec[2] for vec in eigVec1]])

hopfieldCoeff21 = np.real([component*np.conj(component) for component in [vec[0] for vec in eigVec2]])
hopfieldCoeff22 = np.real([component*np.conj(component) for component in [vec[1] for vec in eigVec2]])
hopfieldCoeff23 = np.real([component*np.conj(component) for component in [vec[2] for vec in eigVec2]])

hopfieldCoeff31 = np.real([component*np.conj(component) for component in [vec[0] for vec in eigVec3]])
hopfieldCoeff32 = np.real([component*np.conj(component) for component in [vec[1] for vec in eigVec3]])
hopfieldCoeff33 = np.real([component*np.conj(component) for component in [vec[2] for vec in eigVec3]])

hopfieldCoeff = [[hopfieldCoeff11, hopfieldCoeff12, hopfieldCoeff13],
                 [hopfieldCoeff21, hopfieldCoeff22, hopfieldCoeff23],
                 [hopfieldCoeff31, hopfieldCoeff32, hopfieldCoeff33]]



########################################
##                                    ##
##            Generate plots          ##
##                                    ##
########################################

if not os.path.exists('plots'):
    os.makedirs('plots')


################### Make the eigenvector plots separately ##################################
## 3 for real part and 3 for complex part
# one eigenmode per plot, three components per plot
VecFigs = [[go.Figure(), go.Figure(), go.Figure()], [go.Figure(), go.Figure(), go.Figure()]]
VecFigsNames = [['Eigenmode1_realPart',
                 'Eigenmode2_realPart',
                 'Eigenmode3_realPart'],
                ['Eigenmode1_complexPart',
                 'Eigenmode2_complexPart',
                 'Eigenmode3_complexPart']]
VecFigsTitles = [['Eigenmode 1 (real part)',
                  'Eigenmode 2 (real part)',
                  'Eigenmode 3 (real part)'],
                 ['Eigenmode 1 (complex part)',
                  'Eigenmode 2 (complex part)',
                  'Eigenmode 3 (complex part)']]
componentColors = ['red', 'green', 'blue']
for i in range(len(VecFigs)):
    if (i==0):
        y_data = np.real([eigVec1, eigVec2, eigVec3])
    if (i==1):
        y_data = np.imag([eigVec1, eigVec2, eigVec3])
    for j in range(len(VecFigs[i])):
        fig = VecFigs[i][j]
        fig.add_trace(go.Scatter(x=DelVal, y=[vec[0] for vec in y_data[j]],
                                mode='markers', name='cavity 1'))
        fig.add_trace(go.Scatter(x=DelVal, y=[vec[1] for vec in y_data[j]],
                                mode='markers', name='cavity 2'))
        fig.add_trace(go.Scatter(x=DelVal, y=[vec[2] for vec in y_data[j]],
                                mode='markers', name='exciton'))
        fig.update_layout(title=VecFigsTitles[i][j],
                          xaxis_title='Detuning',
                          showlegend=True)
        fig.write_html('plots/' + VecFigsNames[i][j] + 'VecPlot.html', auto_open=False)

#################### Make 3 eigenvectors in one plot #########################################
## 3x3 components in one plots
# separate real and complex parts
FigNames = ['RealPart', 'ComplexPart']
FigTitles = ['Real part eigenmodes', 'Complex part eigenmodes']
curveNames = ['cavity 1', 'cavity 2', 'exciton']
groupNames = ['eigenmode 1', 'eigenmode 2', 'eigenmode 3']
for i in range(len(VecFigs)):
    if (i==0):
        y_data = np.real([eigVec1, eigVec2, eigVec3])
    if (i==1):
        y_data = np.imag([eigVec1, eigVec2, eigVec3])
    fig = go.Figure()
    for j in range(len(y_data)):
        for k in range(len(y_data)):
            fig.add_trace(go.Scatter(x=DelVal, y=[vec[k] for vec in y_data[j]],
                                    mode='markers', name=curveNames[k],
                                    legendgroup=groupNames[j],
                                    legendgrouptitle_text=groupNames[j],
                                    marker=dict(color=componentColors[k])))
    fig.update_layout(title=FigTitles[i],
                      xaxis_title='Detuning',
                      showlegend=True)
    fig.write_html('plots/' + FigNames[i] + 'VecPlot.html', auto_open=False)

#################### Plot the eigenvalues ##################################################
## make real and complex part separately
# 3 eigenvalues per plot
minRealVertical = min([min(realCompEigVals[0][i]) for i in range(len(realCompEigVals[0]))])
maxRealVertical = max([max(realCompEigVals[0][i]) for i in range(len(realCompEigVals[0]))])
numOfPoints = len(DelVal)
valNames = ['Eigen Val 1', 'Eigen Val 2', 'Eigen Val 3', 'Complex', 'Real']
valColors = ['rgb(15,25,195)', 'rgb(195,5,25)', 'rgb(10,195,25)']
figs = [go.Figure(), go.Figure()]
for j in range(len(['realPart', 'complexPart'])):
    fig = figs[j]
    for i in range(len(['EigenVal1', 'EigenVal2', 'EigenVal3'])):
        fig.add_trace(go.Scatter(x=DelVal, y=realCompEigVals[j][i], mode='lines', name=valNames[i],
                                 line=dict(width=5, color=valColors[i])))
    if j==0:
        fig.add_trace(go.Scatter(x=DelVal, y=df_Eigen['Ecm'].values+(1+r)/2*C,
                                 mode='lines', name="E+rC/sqrt(r)",
                                 line=dict(width=2, color='rgb(150,195,150)', dash='dash')))
        fig.add_trace(go.Scatter(x=DelVal, y=df_Eigen['Ecm'].values-(1+r)/2*C,
                                 mode='lines', name="E-rC/sqrt(r)",
                                 line=dict(width=2, color='rgb(195,150,150)', dash='dash')))
        fig.add_trace(go.Scatter(x=DelVal, y=df_Eigen['XC'], mode='lines', name='XC',
                                 line=dict(width=2, color='rgb(150,150,195)', dash='dash')))
    fig.update_layout(title= valNames[-j-1] + ' part of Eigenvals',
                      xaxis_title='Detuning',
                      yaxis_title= valNames[-j-1] + ' part of Eigenvals',
                      showlegend=True,
                      legend=dict(yanchor="top", y=1, xanchor="left", x=0))
    fig.write_html('plots/' + valNames[-j-1] + 'PartEigen.html', auto_open=True)

############ Plot the Hopfield coefficients ###############################################
## for each eigenmode separately
# 3 components per plot
hopFigNames = ['Cavity 1', 'Cavity 2', 'Exciton']
for j in range(len(hopfieldCoeff)):
    hopFig = go.Figure()
    for i in range(len(hopfieldCoeff[j])):
        hopFig.add_trace(go.Scatter(x=DelVal, y=hopfieldCoeff[j][i],
                                    mode='markers', name=hopFigNames[i]))
    hopFig.update_layout(title= 'Hopfield Coefficients Eigenmode '+ str(j+1),
                         xaxis_title='Detuning',
                         yaxis_title= 'Coeff',
                         showlegend=True)
    hopFig.write_html('plots/' + 'Hopfield_3x3_mode' + str(j+1) + '.html', auto_open=False)

############ Plot all Hopfield coefficients of all eigenmodes in one plot #################
## 3x3 components in one plot
hopFigNames = ['Cavity 1', 'Cavity 2', 'Exciton']
hopColors = ['red', 'green', 'blue']
hopFigAll = go.Figure()
for j in range(len(hopfieldCoeff)):
    for i in range(len(hopfieldCoeff[j])):
        hopFigAll.add_trace(go.Scatter(x=DelVal, y=hopfieldCoeff[j][i],
                                       mode='markers', name=hopFigNames[i],
                                       marker=dict(color=hopColors[i]),
                                       legendgroup='Eigenmode'+str(j+1),
                                       legendgrouptitle_text='Eigenmode'+str(j+1)))
    hopFigAll.update_layout(title= 'Hopfield Coefficients',
                            xaxis_title='Detuning',
                            yaxis_title= 'Coeff',
                            showlegend=True)
hopFigAll.write_html('plots/' + 'Hopfield_3x3_mode_ALL' + '.html', auto_open=False)
