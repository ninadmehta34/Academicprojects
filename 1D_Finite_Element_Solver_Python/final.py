# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 17:36:19 2020

@author: ninad
"""

import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
import math
from scipy.sparse.linalg import cg
#import pdb; pdb.set_trace()

def Shapefunction(n):
    eta = sym.Symbol('eta')
    q = np.linspace(-1, 1, 100)
    n = int(n)
    nodes = n+1
    x = np.linspace(-1, 1, nodes)
    #print(x)
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

def lineElementShapeFunction(a_Eta, a_Degree, a_LocalNode):
    """Compute the shape function value for a line element

    Args:
        a_Eta (float): coordinate in the local coordinate system where shape function is evaluated
        a_Degree (int): the polynomial degree of the finite element/interpolation
        a_LocalNode (int): ID of the local node for which the shape function is evaluated

    Returns:
        basisVal: single float scalar value of the shape function evaluated in local system

    .. _Google Python Style Guide:
    http://google.github.io/styleguide/pyguide.html
    """
    eta = sym.Symbol('eta')
#calling the function
    a, b = Shapefunction(a_Degree)
    #defining the element id at which we are to get values
    index = a_LocalNode - 1
    #getting the shapefunction in terms of eta
    val = a[index]
    der = b[index]
    
    expr2 = val.subs(eta,a_Eta)
    return expr2


def lineElementShapeDerivatives(a_Eta, a_Degree, a_LocalNode):
    """Compute the shape function derivatives for a line element

    Args:
        a_Eta (float): coordinate in the local coordinate system where shape function derivative is evaluated
        a_Degree (int): the polynomial degree of the finite element/interpolation
        a_LocalNode (int): ID of the local node for which the shape function is evaluated

    Returns:
        val: single float scalar value of the shape function derivative evaluated in local system

    .. _Google Python Style Guide:
    http://google.github.io/styleguide/pyguide.html
    """
    eta = sym.Symbol('eta')
#calling the function
    a, b = Shapefunction(a_Degree)
    #defining the element id at which we are to get values
    index = a_LocalNode - 1
    #getting the shapefunction in terms of eta
    der = b[index]
    expr3 = der.subs(eta,a_Eta)
    return expr3
       
# x = np.linspace(-1, 1, n)
# phi = []
#   for i in 
#     for j!=i
#     phi.append(x[])

def lineElementIsoparametricMap(a_ElementNodeCoordinates, a_Degree, a_Eta):
    eta = sym.Symbol('eta')
    a, b = Shapefunction(a_Degree)
    product = []
    sum = 0
    for el1, el2 in zip(a_ElementNodeCoordinates, a):
        product.append(el1*el2)
        #print(product)
    for num in product:
        sum = sum+num
    expr2 = sum.subs(eta,a_Eta)
    return expr2
   

#a,b = lineElementIsoparametricMap([1, 2, 5], 2, 4)
#print(a)
#print(b)

def lineElementMappingGradient(nodal_coordinates, order, coordinate):
    eta = sym.Symbol('eta')
    a, b = Shapefunction(order)
    product = []
    sum = 0
    for el1, el2 in zip(nodal_coordinates, a):
        product.append(el1*el2)
        #print(product)
    for num in product:
        sum = sum+num
    
    derivativeofiso = sum.diff(eta)
    expr3 = derivativeofiso.subs(eta,coordinate)
    return expr3