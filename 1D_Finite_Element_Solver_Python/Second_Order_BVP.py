import sys
#import FEM_1D_Functions
try:
    import numpy as np
    import matplotlib.pyplot as plt
    import sympy as sym
    import math
    from scipy.sparse.linalg import cg
except ImportError:
    sys.exit('Scientific Python Stack: NumPy+SciPy+MatPlotLib Not Installed!')

try:
    from FEM_1D_Functions import *
except ImportError:
    sys.exit('Functions From FEM_1D_Functions Not Successfully Loaded!')


a = 0;
b = 1;
#defining the domain
a_Domain = np.array([a, b]);
#defining the degree, gauss order, and mesh size. these server as the inputs for the entire code.
a_Degree = 1;
a_GaussOrder = 3;
#change the mesh size as desired over here
a_Size = 0.05;
#defining the K value
K = 2.0
#inputting the nodal coordinates
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



a_Func = lambda x: (K*K)*math.cos(3.14*K*x/1.0) + 5.0*(1-K*K)*math.sin(2*3.14*K*x/1.0)
#a_Func = lambda x: (2.0*2.0)*math.cos(3.14*2.0*x/1.0) + 5.0*(1-2.0*2.0)*math.sin(2*3.14*2.0*x/1.0)

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

def applyDirichlet(k_g, f_g, a_DirNodes, a_DirVals):
   
    for i in a_DirNodes:
        
        u0 = a_DirVals[a_DirNodes[i]]
        
        #Account for Dirichlet BC in global loading vector

        f_g = np.add(f_g,-((k_g[:][i])*u0))
        
        
       
        #Set the values in the global stiffness matrix
        k_g[:,i] = 0
        k_g[i,:] = 0
        k_g[i,i] = 1
        f_g[0][i] = u0
       
    return k_g, f_g




