# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 17:38:57 2020

@author: ninad
"""

import numpy as np
import math

from FEM_1D_Functions import generateMeshNodes, generateMeshConnectivity, lineElementShapeFunction, lineElementShapeDerivatives, lineElementIsoparametricMap, lineElementMappingGradient, getGaussQuadratureWeights, getGaussQuadraturePoints


a = 0;
b = 1;
a_Domain = np.array([a, b]);
a_Degree = 1;
a_GaussOrder = 1;
a_Size = 0.1;
K = 2



conn, f, c = generateMeshNodes(a_Domain, a_Degree, a_Size)
xy = generateMeshConnectivity(conn, a_Degree)

def computeElementStiffness(a_ElementNodes, a_GaussOrder, a_Degree):
    order_of_matrix = a_Degree + 1
    K_matrix = np.zeros(shape=(order_of_matrix, order_of_matrix))
    
    for i in range(0, order_of_matrix):
        for j in range(0, order_of_matrix):
            weight = getGaussQuadratureWeights(a_GaussOrder)
            eta_value = getGaussQuadraturePoints(a_GaussOrder)
            for index in range(0, a_GaussOrder):
                derivative1 = lineElementShapeDerivatives(eta_value[index], a_Degree, i)
                #print(derivative1)
                derivative2 = lineElementShapeDerivatives(eta_value[index], a_Degree, j)
                #print(derivative2)
                K_matrix[i][j] = K_matrix[i][j] + ((weight[index] * derivative1 * derivative2)/(lineElementMappingGradient(a_ElementNodes, a_Degree, eta_value[index])))
    return K_matrix

#input your arguements over here 
#c are the element nodal coordinates
tav = computeElementStiffness(c, a_Degree, a_GaussOrder)
print(tav)


def a_Func(x):
     return x
     #return ((K*K)*math.cos(math.pi*K*x) + 5*(1-K*K)*math.sin(2*math.pi*K*x))
     #return ((4 * math.cos(math.pi * 2 * x)) - (15 * math.sin(4 * math.pi * x )))

def computeElementLoading(a_Func, a_ElementNodes, a_GaussOrder, a_Degree):
    #print(a_ElementNodes)
    order_of_matrix = a_Degree + 1
    F = np.zeros(shape=(1, order_of_matrix))
    weight = getGaussQuadratureWeights(a_GaussOrder)
    eta_value = getGaussQuadraturePoints(a_GaussOrder)
    
    for j in range(0, order_of_matrix):
        for index in range(0, a_GaussOrder):
            function_value = a_Func(lineElementIsoparametricMap(a_ElementNodes, a_Degree, eta_value[index]))
            #function_value = lineElementIsoparametricMap(a_ElementNodes, a_Degree, a_func(eta_value[index]))
            #function_value = a_func(a_ElementNodes, a_Degree, eta_value[index])
            phi_value = lineElementShapeFunction(eta_value[index], a_Degree, j)
            map_grad = lineElementMappingGradient(a_ElementNodes, a_Degree, eta_value[index])
            F[0][j] = F[0][j] + (weight[index] * function_value * phi_value * map_grad)
    
    return F

#input your arguements for compute loading in the function below 
#c, here, is for element nodal coordinates
s = computeElementLoading(a_Func, c, a_Degree, a_GaussOrder)
print(s)

