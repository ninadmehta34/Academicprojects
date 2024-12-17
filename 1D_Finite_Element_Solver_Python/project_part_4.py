# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 17:41:58 2020

@author: ninad
"""

import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
import math 


a = 0;
b = 1;
a_Domain = np.array([a, b]);
a_Degree = 1;
a_GaussOrder = 1;
a_Size = 0.5;
K = 2
from FEM_1D_Functions import generateMeshNodes, generateMeshConnectivity, Shapefunction, lineElementShapeFunction, lineElementShapeDerivatives, lineElementIsoparametricMap, lineElementMappingGradient, getGaussQuadratureWeights, getGaussQuadraturePoints
from Second_Order_BVP import computeElementStiffness, computeElementLoading 


conn, f, c = generateMeshNodes(a_Domain, a_Degree, a_Size)
xy = generateMeshConnectivity(conn, a_Degree)


def a_Func(x):
     return x
    

def assembleGlobalStiffness(a_nodes, a_Connectivity, a_GaussOrder, a_Degree):
    n = len(a_nodes)
    dict = {}
    K_matrix = []
    K_g = np.zeros(shape=(n, n))
    m = a_Degree + 1
    i = 0
    while i < len(a_nodes) - 1:
        dict[i] = []
        dict[i].append(a_nodes[i:i+a_Degree+1])
        i = i + a_Degree
    for key, value in dict.items():
        a_ElementNodes = value[0]
        K_matrix.append(computeElementStiffness(a_ElementNodes, a_GaussOrder, a_Degree))
    num = len(a_Connectivity)
    for e in range(num):
        for i in range(m):
            for j in range(m):
                #print(K_matrix[e][i][j])
                #print([a_Connectivity[e][i]])
                #print([a_Connectivity[e][j]])
                K_g[a_Connectivity[e][i]][a_Connectivity[e][j]] =  K_g[a_Connectivity[e][i]][a_Connectivity[e][j]] + K_matrix[e][i][j]
    return K_g

globalstiff = assembleGlobalStiffness(c, xy, a_Degree, a_GaussOrder)
print(globalstiff)

def assembleGlobalLoading(a_Func, a_nodes, a_Connectivity, a_GaussOrder, a_Degree):            
    n = len(a_nodes)
    dict = {}
    F_matrix = []
    F_g = np.zeros(shape=(n,1))

    m = a_Degree + 1
    i = 0
    while i < len(a_nodes) - 1:
        dict[i] = []
        dict[i].append(a_nodes[i:i+a_Degree+1])
        i = i + a_Degree
 
    for key, value in dict.items():
        
        a_ElementNodeCoordinates = value[0]
  
        F_matrix.append(computeElementLoading(a_Func, a_ElementNodeCoordinates, a_GaussOrder, a_Degree)[0])

    num = len(a_Connectivity)
    
    for e in range(num):
        for i in range(m): 
           
            F_g[a_Connectivity[e][i]] = F_g[a_Connectivity[e][i]] + F_matrix[e][i]
         
    return np.transpose(F_g)


#input your arguements in the following for getting the global loading vector
#c = element nodal coordinates while xy is for connectivity matrix
globalload = assembleGlobalLoading(a_Func, [0, 0.5, 1], [[0, 1], [1, 2]], a_GaussOrder, a_Degree)
print('this is global load')
print(globalload)


