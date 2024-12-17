# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 17:00:40 2020

@author: ninad
"""
import numpy as np
import matplotlib.pyplot as plt
import sympy as sym


a = 0;
b = 1;
a_Domain = np.array([a, b]);
a_Degree = 2;
a_GaussOrder = 1;
a_Size = 0.5;
K = 2
from FEM_1D_Functions import generateMeshNodes, generateMeshConnectivity, Shapefunction, lineElementShapeFunction, lineElementShapeDerivatives, lineElementIsoparametricMap, lineElementMappingGradient, plotMesh, plotSolutions, getGaussQuadratureWeights, getGaussQuadraturePoints



def referenceElementNodes(a_Degree):
    n = a_Degree
    n = int(n)
    nodes = n+1
    refnodes = np.linspace(-1, 1, nodes)
    return refnodes
referencenodes = referenceElementNodes(a_Degree)
print(referencenodes)


def Shapefunction(n):
    eta = sym.Symbol('eta')
    q = np.linspace(-1, 1, 100)
    n = int(n)
    nodes = n+1
    x = np.linspace(-1, 1, nodes)
    values = []
    for i in range(nodes):
       for j in range(nodes):
           if j!=i:
            values.append((eta-x[j])/(x[i]-x[j]))    
       
    #print(values)
    #print(shapefunction)
    shape = []
    shape_dr = []

    for i in range(0, len(values), n):
       multiplier = 1
       #print(multiplier)
       for j in range(i, i+n):
           multiplier = multiplier*values[j]
       shape.append(multiplier)
    for sh in shape:
     derivativeshape = sh.diff(eta)
     shape_dr.append(derivativeshape)
     
    #print(shape)    
    #print(shape_dr)
    return(shape, shape_dr)
#insert your inputs here to get the value of shape function and its derivative at a particular node
#Shapefunc,Shapefuncderivative = Shapefunction(a_Degree)
#print(Shapefunc)
#print(Shapefuncderivative)
    


def lineElementShapeFunction(a_Eta, a_Degree, a_LocalNode):

    eta = sym.Symbol('eta')
#calling the function
    a, b = Shapefunction(a_Degree)
    #defining the element id at which we are to get values
    #index = a_LocalNode - 1
    index = a_LocalNode
    #getting the shapefunction in terms of eta
    val = a[index]
    der = b[index]
    
    expr2 = val.subs(eta,a_Eta)
    return expr2

#define the value of a_Eta
#define the value of a_Local node
#shapefunctionvalue = lineElementShapeFunction(a_Eta, a_Degree, a_LocalNode)

def lineElementShapeDerivatives(a_Eta, a_Degree, a_LocalNode):

    eta = sym.Symbol('eta')
    #calling the function
    a, b = Shapefunction(a_Degree)
    #defining the element id at which we are to get values
    index = a_LocalNode 
    #getting the shapefunction in terms of eta
    der = b[index]
    expr3 = der.subs(eta,a_Eta)
    return expr3

#insert yout inputs in the function below to get the values out of it.
#shapefunctionderivative = lineElementShapeDerivatives(a_Eta, a_Degree, a_LocalNode)

def lineElementIsoparametricMap(a_ElementNodeCoordinates, a_Degree, a_Eta):

    eta = sym.Symbol('eta')
    a, b = Shapefunction(a_Degree)
    
    product = []
    sum = 0
    for el1, el2 in zip(a_ElementNodeCoordinates, a):
        
        product.append(el1*el2)
        
    for num in product:
        sum = sum+num
    expr2 = sum.subs(eta,a_Eta)
    return expr2
#insert your inputs in the lineisoparametricvalue to get output. 
#note that the first input is elementnodalcoordinates stored in the variable c

#lineisoparametricvalue = lineElementIsoparametricMap(c, a_Degree, a_Eta)


def lineElementMappingGradient(a_ElementNodeCoordinates, a_Degree, a_Eta):

    eta = sym.Symbol('eta')
    a, b = Shapefunction(a_Degree)
    product = []
    sum = 0
    for el1, el2 in zip(a_ElementNodeCoordinates, a):
        product.append(el1*el2)
        #print(product)
    for num in product:
        sum = sum+num
    
    derivativeofiso = sum.diff(eta)
    expr3 = derivativeofiso.subs(eta,a_Eta)
    return expr3

#insert your inputs to get the value in the lineisoparametricmappinggradientvalue.
##note that the first input is elementnodalcoordinates stored in the variable c
#lineisoparametricmappinggradientvalue  = lineElementMappingGradient(c, a_Degree, a_Eta):
    
#getting the plots for the lineElementShapeFunction

#a_Degree = int(input("Enter the polynomial order: "))
n = np.linspace(-1,1,10)

if a_Degree == 1:

    temp1 = []; temp2 = []; temp3 = []; temp4 = []
    for i in n:
        l1 = lineElementShapeFunction(i, a_Degree, 0)
        l2 = lineElementShapeFunction(i, a_Degree, 1)
        lx1 = lineElementShapeDerivatives(i, a_Degree, 0)
        lx2 =lineElementShapeDerivatives(i, a_Degree, 1)
        temp1.append(l1); temp2.append(l2)
        temp3.append(lx1); temp4.append(lx2)
    plt.subplot(121)
    plt.plot(n,temp1); plt.plot(n,temp2)
    plt.subplot(122)
    plt.plot(n,temp3); plt.plot(n,temp4)
    
    
elif a_Degree == 2:
    n = np.linspace(-1,1,10)
    temp1 = []; temp2 = []; temp3 = []
    temp4 = []; temp5 = []; temp6 = []
    for i in n:
        l1 = lineElementShapeFunction(i, a_Degree, 0)
        l2 = lineElementShapeFunction(i, a_Degree, 1)
        l3 = lineElementShapeFunction(i, a_Degree, 2)
        lx1 = lineElementShapeDerivatives(i, a_Degree, 0)
        lx2 = lineElementShapeDerivatives(i, a_Degree, 1)
        lx3 = lineElementShapeDerivatives(i, a_Degree, 2)
        temp1.append(l1); temp2.append(l2); temp3.append(l3)
        temp4.append(lx1); temp5.append(lx2); temp6.append(lx3)
    plt.subplot(121)    
    plt.plot(n,temp1); plt.plot(n,temp2); plt.plot(n,temp3)
    plt.subplot(122)    
    plt.plot(n,temp4); plt.plot(n,temp5); plt.plot(n,temp6)
    
    
elif a_Degree == 3:
    n = np.linspace(-1,1,10)
    temp1 = []; temp2 = []; temp3 = []; temp4 = []
    temp5 = []; temp6 = []; temp7 = []; temp8 = []
    for i in n:
        l1 = lineElementShapeFunction(i, a_Degree, 0)
        l2 = lineElementShapeFunction(i, a_Degree, 1)
        l3 = lineElementShapeFunction(i, a_Degree, 2)
        l4 = lineElementShapeFunction(i, a_Degree, 3)
        lx1 = lineElementShapeDerivatives(i, a_Degree, 0)
        lx2 = lineElementShapeDerivatives(i, a_Degree, 1)
        lx3 = lineElementShapeDerivatives(i, a_Degree, 2)
        lx4 = lineElementShapeDerivatives(i, a_Degree, 3)
        temp1.append(l1); temp2.append(l2); temp3.append(l3); temp4.append(l4)
        temp5.append(lx1); temp6.append(lx2); temp7.append(lx3); temp8.append(lx4)
    plt.subplot(121)
    plt.plot(n,temp1); plt.plot(n,temp2); plt.plot(n,temp3); plt.plot(n,temp4)
    plt.subplot(122)
    plt.plot(n,temp5); plt.plot(n,temp6); plt.plot(n,temp7); plt.plot(n,temp8)
    
    
elif a_Degree == 4:
    n = np.linspace(-1,1,10)
    temp1 = []; temp2 = []; temp3 = []; temp4 = []; temp5 = []
    temp6 = []; temp7 = []; temp8 = []; temp9 = []; temp10 = []
    for i in n:
        l1 = lineElementShapeFunction(i, a_Degree, 0)
        l2 = lineElementShapeFunction(i, a_Degree, 1)
        l3 = lineElementShapeFunction(i, a_Degree, 2)
        l4 = lineElementShapeFunction(i, a_Degree, 3) 
        l5 = lineElementShapeFunction(i, a_Degree, 4)       
        lx1 = lineElementShapeDerivatives(i, a_Degree, 0)
        lx2 = lineElementShapeDerivatives(i, a_Degree, 1)
        lx3 = lineElementShapeDerivatives(i, a_Degree, 2)
        lx4 = lineElementShapeDerivatives(i, a_Degree, 3)
        lx5 = lineElementShapeDerivatives(i, a_Degree, 4)
        temp1.append(l1); temp2.append(l2); temp3.append(l3); temp4.append(l4); temp5.append(l5)
        temp6.append(lx1); temp7.append(lx2); temp8.append(lx3); temp9.append(lx4); temp10.append(lx5)
    plt.subplot(121)
    plt.plot(n,temp1); plt.plot(n,temp2); plt.plot(n,temp3); plt.plot(n,temp4); plt.plot(n,temp5)
    plt.subplot(122)
    plt.plot(n,temp6); plt.plot(n,temp7); plt.plot(n,temp8); plt.plot(n,temp9); plt.plot(n,temp10)
    

plt.show()
