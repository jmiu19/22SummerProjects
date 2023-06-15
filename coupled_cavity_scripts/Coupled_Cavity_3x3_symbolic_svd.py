import numpy as np
import pandas as pd
import sympy as sym
from sympy import I
import plotly.graph_objects as go
import os

########################################
##                                    ##
##       Setting up the model         ##
##                                    ##
########################################

X, G, R = sym.symbols('X G R')
C = sym.symbols('C')

H = sym.Matrix([[    I*G,   R,   1     ],
                [      R,   X,   R     ],
                [      1,   R,   -I*G  ]])

# (P,D) = H.diagonalize()
k = sym.symbols('k')
charaPoly = sym.det(H-k*sym.eye(3))
eigs = sym.solveset(charaPoly, k)

H.singular_values()


E = sym.Matrix([[      0,   R,   C      ],
                [      R,   X,   R      ],
                [      C,   R,   0      ]])

R = sym.Matrix([[    I*G,   0,   0      ],
                [      0,   0,   0      ],
                [      0,   0,   -I*G   ]])

kE, kR = sym.symbols('kE kR')
charaPolyE = sym.det(E-kE*sym.eye(3))
charaPolyR = sym.det(R-kR*sym.eye(3))
eigsE = sym.solveset(charaPolyE, kE)
eigsR = sym.solveset(charaPolyR, kR)

print(eigsE)
