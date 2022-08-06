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

<<<<<<< HEAD
amax = 800 # number of iteration when computing eigenvalues
=======
amax = 500 # number of iteration when computing eigenvalues
>>>>>>> 501c2673fea82708740fd2746e76e0485560dec2

## Default constants
C = 8 # Coupling
detuning = 7.5
XC = 1650
Rabic = 10
RabiC = 10
G= 0 # Gain
Gu = 30
L='G' # Loss

# # ask user for constants input
# C = float(input('Coupling coefficient between photonic cavities: '))
# Rabic = float(input('Coupling coefficient between photonic cavity 1 and exciton: '))
# RabiC = float(input('Coupling coefficient between photonic cavity 2 and exciton: '))
# detuning = float(input('Detuning of the photonic cavities: '))
# XC = float(input('Energy of the exciton: '))
# G = float(input('Lowest gain: '))
# Gu = float(input('Highest gain: '))
# L = input('Loss of the cavity: (input G if setting loss = gain)   ')

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
                                 'XC'])

df_THT_Eigen = pd.DataFrame(columns=['eigVal',
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

initialG = G
initialL = L

# Compute eigenvals for H
for a in range(0,amax+1):
    G = initialG + (Gu-initialG)/(amax)*(a)
    if (initialL == 'G'):
        L = G
    else:
        L = float(initialL)
    H = np.array([[  Ecm-1j*G,            C,    RabiC/2],
                  [         C,     Ecm+1j*L,    Rabic/2],
                  [   RabiC/2,      Rabic/2,         XC]])
    eigVal, eigVec = np.linalg.eig(H)
    # append eigenvals and eigenvectors to dataframe
    df_Eigen.loc[len(df_Eigen.index)] = [eigVal, eigVec, G, L, C, RabiC/2, Rabic/2, Ecm, XC]

<<<<<<< HEAD
# Compute eigenvals for THT
=======

>>>>>>> 501c2673fea82708740fd2746e76e0485560dec2
for a in range(0,amax+1):
    G = initialG + (Gu-initialG)/(amax)*(a)
    if (initialL == 'G'):
        L = G
    else:
        L = float(initialL)
    THT = np.array([[              Ecm+C,    1j*G,    RabiC/(np.sqrt(2))],
                    [               1j*G,   Ecm-C,                     0],
                    [ RabiC/(np.sqrt(2)),       0,                    XC]])
    eigVal, eigVec = np.linalg.eig(THT)
    # append eigenvals and eigenvectors to dataframe
    df_THT_Eigen.loc[len(df_THT_Eigen.index)] = [eigVal, eigVec, G, L, C, RabiC/2, Rabic/2, Ecm, XC]

##### NOTE #####
# the output eigenvectors is a 3x3 numpy array matrix, each column is an eigenvector,
# but the numpy matrix itself is an array, or a list of rows
# that is eigVec[0] gives the first row, not the first column

## output the dataframe
df_Eigen.to_csv('result.csv')


##############Analytical Eigenvalues###############################
R = Rabic/2
<<<<<<< HEAD
E = Ecm
X = XC
kappa, sigma, zeta, mu1, mu2 = [], [], [], [], []
AnaVal1, AnaVal2, AnaVal3 = [], [], []
diffLMN = []
diff = []

# Compute the eigenvals and parameters anatically
for a in range(0,amax+1):
    G = initialG + (Gu-initialG)/(amax)*(a)
=======
kappa = []
sigma = []
zeta = []
mu1 = []
mu2 = []
AnaVal1 = []
AnaVal2 = []
AnaVal3 = []


for a in range(0,amax+1):
    G = initialG + (Gu-initialG)/(amax)*(a)

    kappa.append(((C*C)-(Ecm*Ecm)-(G*G)+(2*R*R)-(2*Ecm*XC))/(3) + ((2*Ecm+XC)**2)/9)
    sigma.append(C*R*R-Ecm*R*R+((Ecm*Ecm*XC+G*G*XC-C*C*XC)/2)-((2*Ecm+XC)*(-C*C+Ecm*Ecm+2*Ecm*XC+G*G-2*R*R))/(6)+((2*Ecm+XC)**3)/(27))
    zeta.append((sigma[-1] + np.sqrt((sigma[-1]*sigma[-1]) - (kappa[-1]**3) + (0*1j)))**(1/3))

    mu1.append(kappa[-1]/zeta[-1] - zeta[-1])
    mu2.append(kappa[-1]/zeta[-1] + zeta[-1])

    AnaVal1 = (2*Ecm+XC)/3 + (kappa[-1]/zeta[-1] + zeta[-1])
    AnaVal2 = (2*Ecm+XC)/3 - (1/2)*(kappa[-1]/zeta[-1] + zeta[-1]) - (np.sqrt(3)/2)*(kappa[-1]/zeta[-1] - zeta[-1])*1j
    AnaVal3 = (2*Ecm+XC)/3 - (1/2)*(kappa[-1]/zeta[-1] + zeta[-1]) + (np.sqrt(3)/2)*(kappa[-1]/zeta[-1] - zeta[-1])*1j
>>>>>>> 501c2673fea82708740fd2746e76e0485560dec2

    L = (C**2) - (E**2) - (G**2) + (2*R*R) - (2*E*X)
    M = 2*E+X
    N = (E*E*X+G*G*X-C*C*X)/2 + C*R*R - E*R*R
    diffLMN.append(N*N + L*M*N/3 + L*L*M*M/36 + 2*N*M*M*M - L*L*L/27 - L*L*M*M/27)

    kappa.append(+ ((C*C)-(Ecm*Ecm)-(G*G)+(2*R*R)-(2*Ecm*XC))/(3)
                 + ((2*Ecm+XC)**2)/9)
    sigma.append(+ C*R*R
                 - Ecm*R*R
                 + ((Ecm*Ecm*XC+G*G*XC-C*C*XC)/2)
                 - ((2*Ecm+XC)*(-C*C+Ecm*Ecm+2*Ecm*XC+G*G-2*R*R))/(6)
                 + ((2*Ecm+XC)**3)/(27))
    zeta.append((+ sigma[-1]
                 + np.sqrt((sigma[-1]*sigma[-1])
                 - (kappa[-1]**3) + (0*1j)))**(1/3))

    mu1.append(kappa[-1]/zeta[-1] - zeta[-1])
    mu2.append(kappa[-1]/zeta[-1] + zeta[-1])

    AnaVal1.append((2*Ecm+XC)/3 + (kappa[-1]/zeta[-1] + zeta[-1]))
    AnaVal2.append((2*Ecm+XC)/3 - (1/2)*(kappa[-1]/zeta[-1] + zeta[-1])
                                - (np.sqrt(3)/2)*(kappa[-1]/zeta[-1] - zeta[-1])*1j)
    AnaVal3.append((2*Ecm+XC)/3 - (1/2)*(kappa[-1]/zeta[-1] + zeta[-1])
                                + (np.sqrt(3)/2)*(kappa[-1]/zeta[-1] - zeta[-1])*1j)

    diff.append(kappa[-1]**3 - sigma[-1]**2)

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


########### THT eigenstuff #############################
THTrealEigVal = np.real(df_THT_Eigen['eigVal'].values.tolist())
THTcompEigVal = np.imag(df_THT_Eigen['eigVal'].values.tolist())
THTgainVal = np.real(df_THT_Eigen['Gain'].values.tolist())

THTrealEigVal1 = [vals[0] for vals in THTrealEigVal]
THTrealEigVal2 = [vals[1] for vals in THTrealEigVal]
THTrealEigVal3 = [vals[2] for vals in THTrealEigVal]

THTcompEigVal1 = [vals[0] for vals in THTcompEigVal]
THTcompEigVal2 = [vals[1] for vals in THTcompEigVal]
THTcompEigVal3 = [vals[2] for vals in THTcompEigVal]

THTrealCompEigVals = [[THTrealEigVal1, THTrealEigVal2, THTrealEigVal3], [THTcompEigVal1, THTcompEigVal2, THTcompEigVal3]]




RMatchV1V2Count, RMatchV1V3Count, RMatchV2V3Count = 0, 0, 0
CMatchV1V2Count, CMatchV1V3Count, CMatchV2V3Count = 0, 0, 0

RMatchV1V2Index, RMatchV1V3Index, RMatchV2V3Index = [], [], []
CMatchV1V2Index, CMatchV1V3Index, CMatchV2V3Index = [], [], []

MatchCount = [[RMatchV1V2Count, RMatchV1V3Count, RMatchV2V3Count],
              [CMatchV1V2Count, CMatchV1V3Count, CMatchV2V3Count]]
MatchIndex = [[RMatchV1V2Index, RMatchV1V3Index, RMatchV2V3Index],
              [CMatchV1V2Index, CMatchV1V3Index, CMatchV2V3Index]]


########### THT eigenstuff #############################
THTrealEigVal = np.real(df_THT_Eigen['eigVal'].values.tolist())
THTcompEigVal = np.imag(df_THT_Eigen['eigVal'].values.tolist())
THTgainVal = np.real(df_THT_Eigen['Gain'].values.tolist())

THTrealEigVal1 = [vals[0] for vals in THTrealEigVal]
THTrealEigVal2 = [vals[1] for vals in THTrealEigVal]
THTrealEigVal3 = [vals[2] for vals in THTrealEigVal]

THTcompEigVal1 = [vals[0] for vals in THTcompEigVal]
THTcompEigVal2 = [vals[1] for vals in THTcompEigVal]
THTcompEigVal3 = [vals[2] for vals in THTcompEigVal]

THTrealCompEigVals = [[THTrealEigVal1, THTrealEigVal2, THTrealEigVal3], [THTcompEigVal1, THTcompEigVal2, THTcompEigVal3]]


## Extract the eigenvectors
eigVecsSet = df_Eigen['eigVec']
gainVal = np.real(df_Eigen['Gain'].values.tolist())

##### NOTE ##### refers to notes in the Compute Eigenvals section
# the output eigenvectors is a 3x3 numpy array matrix, each column is an eigenvector,
# but the numpy matrix itself is an array, or a list of rows
# that is eigVec[0] gives the first row, not the first column

eigVec1 = [[vecs[0][0], vecs[1][0], vecs[2][0]] for vecs in eigVecsSet]
eigVec2 = [[vecs[0][1], vecs[1][1], vecs[2][1]] for vecs in eigVecsSet]
eigVec3 = [[vecs[0][2], vecs[1][2], vecs[2][2]] for vecs in eigVecsSet]

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

orderParameter2 = []
# the second way of compute order parameter
for i in range(len(realCompEigVals[1][1])):
    orderParameter2.append(abs(realCompEigVals[1][0][i])+abs(realCompEigVals[1][1][i])+abs(realCompEigVals[1][2][i]))

# the first way of computing order parameter
orderParameter1 = []
eigenVecs = [eigVec1, eigVec2, eigVec3]
for i in range(len(eigVec1)):
    eta = 0
    for j in range(len(eigenVecs)):
        PTvec = np.transpose(np.matmul(P, np.conj([[eigenVecs[j][i][0]],
                                                   [eigenVecs[j][i][1]],
                                                   [eigenVecs[j][i][2]]])))[0]
        eta = eta + abs(+ np.dot(np.conj(eigenVecs[j][i]), eigenVecs[j][i])
                        - np.dot(np.conj(PTvec), eigenVecs[j][i]))
    orderParameter1.append(eta)


########################################
##                                    ##
##       Find exceptional points      ##
##                                    ##
########################################

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
            excepGain = initialG + (Gu-initialG)/(amax)*MatchIndex[0][i][k]
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

if not os.path.exists('plots'):
    os.makedirs('plots')

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
        fig.add_trace(go.Scatter(x=gainVal, y=[vec[0] for vec in y_data[j]],
                                mode='markers', name='cavity 1'))
        fig.add_trace(go.Scatter(x=gainVal, y=[vec[1] for vec in y_data[j]],
                                mode='markers', name='cavity 2'))
        fig.add_trace(go.Scatter(x=gainVal, y=[vec[2] for vec in y_data[j]],
                                mode='markers', name='exciton'))
        fig.update_layout(title=VecFigsTitles[i][j],
                          xaxis_title='Gain of cavity 1',
                          showlegend=True)
        fig.write_html('plots/' + VecFigsNames[i][j] + 'VecPlot.html', auto_open=False)

#################### Make 3 eigenmodes in one plot #########################################
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
            fig.add_trace(go.Scatter(x=gainVal, y=[vec[k] for vec in y_data[j]],
                                    mode='markers', name=curveNames[k],
                                    legendgroup=groupNames[j],
                                    legendgrouptitle_text=groupNames[j],
                                    marker=dict(color=componentColors[k])))
    fig.update_layout(title=FigTitles[i],
                      xaxis_title='Gain of cavity 1',
                      showlegend=True)
    fig.write_html('plots/' + FigNames[i] + 'VecPlot.html', auto_open=False)

#################### Plot the eigenvalues ##################################################
## make real and complex part separately
minRealVertical = min([min(realCompEigVals[0][i]) for i in range(len(realCompEigVals[0]))])
maxRealVertical = max([max(realCompEigVals[0][i]) for i in range(len(realCompEigVals[0]))])
numOfPoints = len(gainVal)
valNames = ['Eigen Val 1', 'Eigen Val 2', 'Eigen Val 3', 'Complex', 'Real']
valColors = ['rgb(15,25,195)', 'rgb(195,5,25)', 'rgb(10,195,25)']
figs = [go.Figure(), go.Figure()]
for j in range(len(['realPart', 'complexPart'])):
    fig = figs[j]
    for i in range(len(['EigenVal1', 'EigenVal2', 'EigenVal3'])):
        fig.add_trace(go.Scatter(x=gainVal, y=realCompEigVals[j][i], mode='markers', name=valNames[i],
                                 marker=dict(size=10, color=valColors[i])))
        fig.add_trace(go.Scatter(x=[initialG + (Gu-initialG)/(amax)*(val-1) for val in excepIndex[i]],
                                 y=[realCompEigVals[j][i][t] for t in excepIndex[i]],
                                 mode='markers', name='Exceptional Point', marker=dict(size=19, color='black'),
                                 legendgroup='Exceptional Points', legendgrouptitle_text='Exceptional Points',
                                 legendrank=1001))
    if (j==0):
        fig.add_trace(go.Scatter(x=gainVal, y=np.linspace(Ecm, Ecm, num=numOfPoints),
                                 mode='lines', name='Ecm', legendgroup='Parameters',
                                 legendgrouptitle_text='Parameters',
                                 line=dict(color='rgb(255,175,175)', width=2, dash='dot')))
        fig.add_trace(go.Scatter(x=gainVal, y=np.linspace(XC, XC, num=numOfPoints),
                                 mode='lines', name='Xc', legendgroup='Parameters',
                                 legendgrouptitle_text='Parameters',
                                 line=dict(color='rgb(200,200,255)', width=2, dash='dot')))
        fig.add_trace(go.Scatter(x=gainVal, y=np.linspace(Ecm+C, Ecm+C, num=numOfPoints),
                                 mode='lines', name='Ecm + C', legendgroup='Parameters',
                                 legendgrouptitle_text='Parameters',
                                 line=dict(color='rgb(150,255,185)', width=2, dash='dot')))
        fig.add_trace(go.Scatter(x=gainVal, y=np.linspace(Ecm-C, Ecm-C, num=numOfPoints),
                                 mode='lines', name='Ecm - C', legendgroup='Parameters',
                                 legendgrouptitle_text='Parameters',
                                 line=dict(color='rgb(150,255,185)', width=2, dash='dot')))
        fig.add_trace(go.Scatter(x=[C for index in range(len(gainVal))],
                                 y=np.linspace(minRealVertical, maxRealVertical, num=numOfPoints),
                                 mode='lines', name='C', legendgroup='Parameters',
                                 legendgrouptitle_text='Parameters',
                                 line=dict(color='rgb(240,240,170)', width=2, dash='dot')))
    fig.update_layout(title= valNames[-j-1] + ' part of Eigenvals',
                      xaxis_title='Gain of Band 1',
                      yaxis_title= valNames[-j-1] + ' part of Eigenvals',
                      showlegend=True)
    fig.write_html('plots/' + valNames[-j-1] + 'PartEigen.html', auto_open=True)

############ Plot the Hopfield coefficients ###############################################
## for each eigenmode separately
hopFigNames = ['Gainy Cavity', 'Lossy Cavity', 'Exciton']
for j in range(len(hopfieldCoeff)):
    hopFig = go.Figure()
    for i in range(len(hopfieldCoeff[j])):
        hopFig.add_trace(go.Scatter(x=gainVal, y=hopfieldCoeff[j][i],
                                    mode='markers', name=hopFigNames[i]))
    hopFig.update_layout(title= 'Hopfield Coefficients Eigenmode '+ str(j+1),
                         xaxis_title='Gain of Band 1',
                         yaxis_title= 'Coeff',
                         showlegend=True)
    hopFig.write_html('plots/' + 'Hopfield_3x3_mode' + str(j+1) + '.html', auto_open=False)

############ Plot all Hopfield coefficients of all eigenmodes in one plot #################
hopFigNames = ['Gainy Cavity', 'Lossy Cavity', 'Exciton']
hopColors = ['red', 'green', 'blue']
hopFigAll = go.Figure()
for j in range(len(hopfieldCoeff)):
    for i in range(len(hopfieldCoeff[j])):
        hopFigAll.add_trace(go.Scatter(x=gainVal, y=hopfieldCoeff[j][i],
                                       mode='markers', name=hopFigNames[i],
                                       marker=dict(color=hopColors[i]),
                                       legendgroup='Eigenmode'+str(j+1),
                                       legendgrouptitle_text='Eigenmode'+str(j+1)))
    hopFigAll.update_layout(title= 'Hopfield Coefficients',
                            xaxis_title='Gain of Band 1',
                            yaxis_title= 'Coeff',
                            showlegend=True)
hopFigAll.write_html('plots/' + 'Hopfield_3x3_mode_ALL' + '.html', auto_open=True)




############# Make the phase transition plot ##############################################
phaseTransPlot = go.Figure()
phaseTransPlot.add_trace(go.Scatter(x=gainVal, y=orderParameter1,
                                    mode='markers', name='1st method'))
phaseTransPlot.add_trace(go.Scatter(x=gainVal, y=orderParameter2,
                                    mode='markers', name='2nd method'))
phaseTransPlot.update_layout(title= 'Order Paramter',
                             xaxis_title='Gain of Band 1',
                             yaxis_title= 'eta',
                             showlegend=True)
phaseTransPlot.write_html('plots/' + 'orderParameter.html', auto_open=False)


<<<<<<< HEAD
############## Plot analytic solution parameters ##########################################
AnaParaPlot = go.Figure()
AnaParaPlot.add_trace(go.Scatter(x=gainVal, y=np.real(diff),
                                 mode='markers', name='kappa3 - sigma2'))
AnaParaPlot.add_trace(go.Scatter(x=gainVal, y=np.imag(diffLMN),
                                 mode='markers', name='kappa3 - sigma2 (LMN)'))
# AnaParaPlot.add_trace(go.Scatter(x=gainVal, y=np.real(mu2),
#                                  mode='markers', name='real of mu2 (kappa/zeta + zeta)'))
# AnaParaPlot.add_trace(go.Scatter(x=gainVal, y=np.imag(mu2),
#                                  mode='markers', name='imag of mu2'))
# AnaParaPlot.add_trace(go.Scatter(x=gainVal, y=kappa,
#                                  mode='markers', name='kappa'))
AnaParaPlot.add_trace(go.Scatter(x=gainVal, y=[np.real(kappa[i]-(zeta[i]**2)) for i in range(len(kappa))],
                                 mode='markers', name='kappa - zeta2 real'))
AnaParaPlot.add_trace(go.Scatter(x=gainVal, y=[np.imag(kappa[i]-zeta[i]) for i in range(len(kappa))],
                                 mode='markers', name='kappa - zeta2 image'))
# AnaParaPlot.add_trace(go.Scatter(x=gainVal, y=np.imag(zeta),
#                                  mode='markers', name='imag of zeta'))
AnaParaPlot.update_layout(title= 'Analytic solution parameters',
                          xaxis_title='Gain of Band 1',
                          yaxis_title= 'vals',
                          showlegend=True)
AnaParaPlot.write_html('plots/' + 'AnalyticSolutionParameters.html', auto_open=False)


############ Plot analytic solution eigenvalues ###############################################
AnaValRPlot = go.Figure()
AnaValRPlot.add_trace(go.Scatter(x=gainVal, y=np.real(AnaVal1),
                                 mode='markers', name='Eigval 1'))
AnaValRPlot.add_trace(go.Scatter(x=gainVal, y=np.real(AnaVal2),
                                 mode='markers', name='Eigval 2'))
AnaValRPlot.add_trace(go.Scatter(x=gainVal, y=np.real(AnaVal3),
                                 mode='markers', name='Eigval 3'))
AnaValRPlot.update_layout(title= 'Analytic solution eigenvalues real part',
                          xaxis_title='Gain of Band 1',
                          yaxis_title= 'vals',
                          showlegend=True)
AnaValRPlot.write_html('plots/' + 'AnalyticRealEigenVals.html', auto_open=False)

AnaValCPlot = go.Figure()
AnaValCPlot.add_trace(go.Scatter(x=gainVal, y=np.imag(AnaVal1),
                                 mode='markers', name='Eigval 1'))
AnaValCPlot.add_trace(go.Scatter(x=gainVal, y=np.imag(AnaVal2),
                                 mode='markers', name='Eigval 2'))
AnaValCPlot.add_trace(go.Scatter(x=gainVal, y=np.imag(AnaVal3),
                                 mode='markers', name='Eigval 3'))
AnaValCPlot.update_layout(title= 'Analytic solution eigenvalues imaginary part',
                          xaxis_title='Gain of Band 1',
                          yaxis_title= 'vals',
                          showlegend=True)
AnaValCPlot.write_html('plots/' + 'AnalyticImaginaryEigenVals.html', auto_open=False)


################# Plot the THT eigenvalues ###################################################
=======
###########################################################################################
## plot analytic solution of eigenvals
phaseTransPlot = go.Figure()
phaseTransPlot.add_trace(go.Scatter(x=gainVal, y=np.real(mu1),
                                    mode='markers', name='real of mu1'))
phaseTransPlot.add_trace(go.Scatter(x=gainVal, y=np.imag(mu1),
                                    mode='markers', name='imag of mu1'))
phaseTransPlot.add_trace(go.Scatter(x=gainVal, y=np.real(mu2),
                                    mode='markers', name='real of mu2'))
phaseTransPlot.add_trace(go.Scatter(x=gainVal, y=np.imag(mu2),
                                    mode='markers', name='imag of mu2'))
phaseTransPlot.update_layout(title= 'Analytic solution parameters',
                             xaxis_title='Gain of Band 1',
                             yaxis_title= 'vals',
                             showlegend=True)
phaseTransPlot.write_html('plots/' + 'AnalyticSolutionParameters.html', auto_open=True)



###########################################################################################
## plot the eigenvalues
>>>>>>> 501c2673fea82708740fd2746e76e0485560dec2
# make real and complex part separately
numOfPoints = len(gainVal)
valNames = ['Eigen Val 1', 'Eigen Val 2', 'Eigen Val 3', 'Complex', 'Real']
valColors = ['rgb(15,25,195)', 'rgb(195,5,25)', 'rgb(10,195,25)']
figs = [go.Figure(), go.Figure()]
<<<<<<< HEAD
=======

>>>>>>> 501c2673fea82708740fd2746e76e0485560dec2
for j in range(len(['realPart', 'complexPart'])):
    fig = figs[j]
    for i in range(len(['EigenVal1', 'EigenVal2', 'EigenVal3'])):
        fig.add_trace(go.Scatter(x=gainVal, y=THTrealCompEigVals[j][i], mode='markers', name=valNames[i],
                                 marker=dict(size=10, color=valColors[i])))
    fig.update_layout(title= valNames[-j-1] + ' part of THT Eigenvals',
                      xaxis_title='Gain of Band 1',
                      yaxis_title= valNames[-j-1] + ' part of THT Eigenvals',
                      showlegend=True)
    fig.write_html('plots/' + valNames[-j-1] + 'PartTHTEigen.html', auto_open=False)
