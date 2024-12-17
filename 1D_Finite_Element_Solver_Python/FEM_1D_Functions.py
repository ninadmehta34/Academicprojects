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
a_GaussOrder = 3;
#change the mesh size as desired over here
a_Size = 0.05;
#defining the K value
K = 2.0

def generateMeshNodes(a_Domain, a_Degree, a_Size, a_ReturnNumElements=True, a_ReturnNumNodes=False):

    domain = np.array([a_Domain[0], a_Domain[1]])
    numElements = (domain[1] - domain[0])/a_Size
    n = 1+a_Degree
    numNodes = (numElements*a_Degree)+1
    numNodes = int(numNodes)
    nodes = np.linspace(domain[0],domain[1],numNodes)
    return [numElements, numNodes, nodes]


#the output of this functions are stored in the following variables:
# conn gives the number of elements, f gives the number of nodes and c gives the element nodal coordinates needed.
#conn, f, c = generateMeshNodes(a_Domain, a_Degree, a_Size)



def generateMeshConnectivity(a_NumElements, a_Degree):

    a_NumElements = int(a_NumElements)
    connectivity = []
    #a_NumElements = 10
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

#xy = generateMeshConnectivity(conn, a_Degree)
#print(xy)


def referenceElementNodes(a_Degree):
    n = a_Degree
    n = int(n)
    nodes = n+1
    refnodes = np.linspace(-1, 1, nodes)
    return refnodes
#referencenodes = referenceElementNodes(a_Degree)



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


def lineElementShapeDerivatives(a_Eta, a_Degree, a_LocalNode):

    eta = sym.Symbol('eta')
    #calling the function
    a, b = Shapefunction(a_Degree)
    #defining the element id at which we are to get values
    #index = a_LocalNode - 1
    index = a_LocalNode 
    #getting the shapefunction in terms of eta
    der = b[index]
    expr3 = der.subs(eta,a_Eta)
    return expr3



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

def plotMesh(a_Nodes, a_Connectivity, a_XLabel=None, a_YLabel=None, a_SaveFigure=None):


    if a_Connectivity.shape[1] == 2:

        plt.plot(a_Nodes, np.zeros_like(a_Nodes), 'rs-')

    else:
        edgeNodes   = np.unique(np.concatenate((a_Nodes[a_Connectivity[:,0]], a_Nodes[a_Connectivity[:,-1]])))
        inNodes     = a_Nodes[np.array([x not in edgeNodes for x in a_Nodes])]
        print(inNodes)
        plt.plot(edgeNodes, np.zeros_like(edgeNodes), 'rs-')
        plt.plot(inNodes, np.zeros_like(inNodes), 'bv')

    plt.show()



def plotSolutions(a_Nodes, a_Connectivity, a_Solution, a_YLabel=None, a_XLabel=None, a_SaveFigure=None):


    if a_Connectivity.shape[1] == 2:

        plt.plot(a_Nodes, np.zeros_like(a_Nodes), 'rs-')

    else:

        edgeNodes   = np.unique(np.concatenate((a_Nodes[a_Connectivity[:,0]], a_Nodes[a_Connectivity[:,-1]])))
        inNodes     = a_Nodes[np.array([x not in edgeNodes for x in a_Nodes])]
        plt.plot(edgeNodes, np.zeros_like(edgeNodes), 'rs-')
        plt.plot(inNodes, np.zeros_like(inNodes), 'bv')

    plt.plot(a_Nodes, a_Solution, 'ms-')

    if a_YLabel is not None:
        plt.ylabel(a_YLabel, fontweight='bold')

    if a_XLabel is not None:
        plt.xlabel(a_XLabel, fontweight='bold')

    if a_SaveFigure is not None:
        plt.savefig(a_SaveFigure, bbox_inches='tight')
        plt.close()
    else:
        plt.show()


#-----------------------------------------------------------------------
# Tabulation of weights and quadrature points for 1D Gaussian Quadrature
#-----------------------------------------------------------------------
def getGaussQuadratureWeights(a_Degree):


    if a_Degree == 1:
        return [2.0]
    elif a_Degree == 2:
        return [1.0, 1.0]
    elif a_Degree == 3:
        return [5.0/9.0, 8.0/9.0, 5.0/9.0]
    elif a_Degree == 4:
        return [0.652145, 0.347855, 0.652145, 0.347855]
    else:
        sys.exit("Only upto 4th degree Gauss quadrature is implemented")


def getGaussQuadraturePoints(a_Degree):


    if a_Degree == 1:
        return [0.0]
    elif a_Degree == 2:
        return [-0.57735, 0.57735]
    elif a_Degree == 3:
        return [-0.774597, 0.0, 0.774597]
    elif a_Degree == 4:
        return [-0.339981, -0.861136, 0.339981, 0.861136]
    else:
        sys.exit("Only upto 4th degree Gauss quadrature is implemented")


