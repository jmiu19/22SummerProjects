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

# ## Default constants
# C = 8 # Coupling
# detuning = 7.5
# XC = 1650
# Rabic = 10
# RabiC = 10
# G= 0 # Gain
# Gu = 30
# L='G' # Loss

# ask user for constants input
C = float(input('value of C: '))
Rabic = float(input('value of 2*Omega: '))
RabiC = Rabic
XC = float(input('energy of the exciton: '))
G = float(input('gain: '))
detuning = float(input('lowest detuning: '))
h_detuning = float(input('highest detuning: '))
L = G

Ecm = XC + detuning

## Create empty dataframe for storing data
df_Eigen = pd.DataFrame(columns=['eigVal',
                                 'eigVec',
                                 'Gain',
                                 'Loss',
                                 'Coupling bt. bands',
                                 'Coupling bt. band1 and exciton',
                                 'Coupling bt. band2 and exciton',
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

# Compute eigenvals for H while sweeping Gain
for a in range(0,amax+1):
    detuning = initialDel + (h_detuning-initialDel)/(amax)*(a)
    Ecm = XC + detuning
    H = np.array([[  Ecm-1j*G,            C,    RabiC/2],
                  [         C,     Ecm+1j*L,    Rabic/2],
                  [   RabiC/2,      Rabic/2,         XC]])
    eigVal, eigVec = np.linalg.eig(H)
    # append eigenvals and eigenvectors to dataframe
    df_Eigen.loc[len(df_Eigen.index)] = [eigVal, eigVec, G, L, C, RabiC/2, Rabic/2, Ecm, XC, detuning]

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

realEigVal1 = [vals[0] for vals in realEigVal]
realEigVal2 = [vals[1] for vals in realEigVal]
realEigVal3 = [vals[2] for vals in realEigVal]

compEigVal1 = [vals[0] for vals in compEigVal]
compEigVal2 = [vals[1] for vals in compEigVal]
compEigVal3 = [vals[2] for vals in compEigVal]

realCompEigVals = [[realEigVal1, realEigVal2, realEigVal3], [compEigVal1, compEigVal2, compEigVal3]]



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
##       Find exceptional points      ##
##                                    ##
########################################
RMatchV1V2Count, RMatchV1V3Count, RMatchV2V3Count = 0, 0, 0
CMatchV1V2Count, CMatchV1V3Count, CMatchV2V3Count = 0, 0, 0

RMatchV1V2Index, RMatchV1V3Index, RMatchV2V3Index = [], [], []
CMatchV1V2Index, CMatchV1V3Index, CMatchV2V3Index = [], [], []

MatchCount = [[RMatchV1V2Count, RMatchV1V3Count, RMatchV2V3Count],
              [CMatchV1V2Count, CMatchV1V3Count, CMatchV2V3Count]]
MatchIndex = [[RMatchV1V2Index, RMatchV1V3Index, RMatchV2V3Index],
              [CMatchV1V2Index, CMatchV1V3Index, CMatchV2V3Index]]

## since the results are computed numerically
## two identical eigenvalues might have a very small numerical difference
## we define the maximal difference two eigenvals can have to be considered as being equal

############ use max difference to define threshold #################
# threshold = 0                                                     #
# for i in range(len(realCompEigVals)):                             #
#     for j in range(len(realCompEigVals[i])):                      #
#         if threshold < max(np.diff(realCompEigVals[i][j]))/3:     #
#             threshold = max(np.diff(realCompEigVals[i][j]))/3     #
#####################################################################
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
            excepDel = initialDel + (h_detuning-initialDel)/(amax)*MatchIndex[0][i][k]
            print('Found potentional exceptional point between '
                  + str(checkIndex[i])
                  + ' at around detuning = '
                  + str(excepDel))
            excepIndex[i].append(MatchIndex[0][i][k])


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
                          xaxis_title='Gain of cavity 1',
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
        fig.add_trace(go.Scatter(x=DelVal, y=realCompEigVals[j][i], mode='markers', name=valNames[i],
                                 marker=dict(size=10, color=valColors[i])))
        fig.add_trace(go.Scatter(x=[initialDel + (h_detuning-initialDel)/(amax)*(val-1) for val in excepIndex[i]],
                                 y=[realCompEigVals[j][i][t] for t in excepIndex[i]],
                                 mode='markers', name='Exceptional Point', marker=dict(size=19, color='black'),
                                 legendgroup='Exceptional Points', legendgrouptitle_text='Exceptional Points',
                                 legendrank=1001))
    fig.update_layout(title= valNames[-j-1] + ' part of Eigenvals',
                      xaxis_title='Detuning',
                      yaxis_title= valNames[-j-1] + ' part of Eigenvals',
                      showlegend=True)
    fig.write_html('plots/' + valNames[-j-1] + 'PartEigen.html', auto_open=True)


#################### Log plot the eigenvalues ##################################################
## make real and complex part separately
# 3 eigenvalues per plot
# log scale the y-axis (real part of the eigenvalues)
minRealVertical = min([min(realCompEigVals[0][i]) for i in range(len(realCompEigVals[0]))])
maxRealVertical = max([max(realCompEigVals[0][i]) for i in range(len(realCompEigVals[0]))])
numOfPoints = len(DelVal)
valNames = ['Eigen Val 1', 'Eigen Val 2', 'Eigen Val 3', 'Complex', 'Log scale Real']
valColors = ['rgb(15,25,195)', 'rgb(195,5,25)', 'rgb(10,195,25)']
figs = [go.Figure(), go.Figure()]
for j in range(len(['realPart', 'complexPart'])):
    fig = figs[j]
    for i in range(len(['EigenVal1', 'EigenVal2', 'EigenVal3'])):
        y_vals = np.array(realCompEigVals[j][i])
        x_vals = np.array(DelVal)
        eps_x = np.array([initialDel + (h_detuning-initialDel)/(amax)*(val-1) for val in excepIndex[i]])
        eps_y = np.array([realCompEigVals[j][i][t] for t in excepIndex[i]])
        if j==0:
            x_vals = np.log(x_vals)
            y_vals = np.log(y_vals)
            eps_x = np.log(eps_x)
            eps_y = np.log([realCompEigVals[j][i][t] for t in excepIndex[i]])
        fig.add_trace(go.Scatter(x=x_vals, y=y_vals,
                                 mode='markers', name=valNames[i],
                                 marker=dict(size=10, color=valColors[i])))
        fig.add_trace(go.Scatter(x=eps_x, y=eps_y,
                                 mode='markers', name='Exceptional Point', marker=dict(size=19, color='black'),
                                 legendgroup='Exceptional Points', legendgrouptitle_text='Exceptional Points',
                                 legendrank=1001))
    fig.update_layout(title= valNames[-j-1] + ' part of Eigenvals',
                      xaxis_title='Detuning',
                      yaxis_title= valNames[-j-1] + ' part of Eigenvals',
                      showlegend=True)
    fig.write_html('plots/' + valNames[-j-1] + 'PartEigen.html', auto_open=True)




############ Plot the Hopfield coefficients ###############################################
## for each eigenmode separately
# 3 components per plot
hopFigNames = ['Gainy Cavity', 'Lossy Cavity', 'Exciton']
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
hopFigNames = ['Gainy Cavity', 'Lossy Cavity', 'Exciton']
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
