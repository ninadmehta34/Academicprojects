# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 16:55:41 2020

@author: ninad
"""

import sys
try:
    import numpy as np
    import matplotlib.pyplot as plt
    import sympy as sym
    import math
    from scipy.sparse.linalg import cg
except ImportError:
    sys.exit('Scientific Python Stack: NumPy+SciPy+MatPlotLib Not Installed!')
#defining the upper and lower limits
a = 0;
b = 1;
#defining the domain
a_Domain = np.array([a, b]);
#defining the degree, gauss order, and mesh size. these server as the inputs for the entire code.
a_Degree = 1;
a_GaussOrder = 1;
a_Size = 0.5;
#defining the K value
K = 2

def generateMeshNodes(a_Domain, a_Degree, a_Size, a_ReturnNumElements=True, a_ReturnNumNodes=False):

    domain = np.array([a_Domain[0], a_Domain[1]])
    print(domain)
    numElements = (domain[1] - domain[0])/a_Size
    n = 1+a_Degree
    numNodes = (numElements*a_Degree)+1
    numNodes = int(numNodes)
    nodes = np.linspace(domain[0],domain[1],numNodes)
    return [numElements, numNodes, nodes]

#the output of this functions are stored in the following variables:
# conn gives the number of elements, f gives the number of nodes and c gives the element nodal coordinates needed.
conn, f, c = generateMeshNodes(a_Domain, a_Degree, a_Size)
print(conn)
print(f)
print(c)


def generateMeshConnectivity(a_NumElements, a_Degree):

    a_NumElements = int(a_NumElements)
    connectivity = []
    element_type = a_Degree + 1
    i = 0
    while len(connectivity) < a_NumElements:
        temp = []
        while len(temp)<element_type:
            temp.append(i)
            if len(temp)!=element_type:
                i = i+1
        connectivity.append(temp)
    return connectivity

xy = generateMeshConnectivity(conn, a_Degree)
print(xy)
