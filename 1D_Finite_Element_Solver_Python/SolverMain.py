import sys
try:
    import numpy as np
    import matplotlib.pyplot as plt
    import math
    from scipy.sparse.linalg import cg
except ImportError:
    sys.exit('Scientific Python Stack: NumPy+SciPy+MatPlotLib Not Installed!')

try:
    from FEM_1D_Functions import *
except ImportError:
    sys.exit('Functions From FEM_1D_Functions Not Successfully Loaded!')

try:
    from Second_Order_BVP import *
except ImportError:
    sys.exit('Functions From Second_Order_BVP Not Successfully Loaded!')

#------------------------
# configure solver inputs
#------------------------
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

#-----------------------------------------
# define the functionS for f
# and the exact solution you have computed
#------------------------------------------
a_Func = lambda x: (K*K)*math.cos(3.14*K*x/1.0) + (5.0*(1-K*K)*math.sin(2*3.14*K*x/1.0))

#exact = lambda x: (-1/(3.14*3.14))*math.cos(K*3.14*x)+(15/(16*3.14*3.14))*math.sin(4*3.14*x)+1.10132*x+0.10132
exact = lambda x: (-1/(math.pi*math.pi))*math.cos(K*math.pi*x)-5*(1-K*K)/(4*math.pi*math.pi*K*K)*math.sin(2*K*math.pi*x)+(1+1/(math.pi*math.pi)*math.cos(math.pi*K))*x+1/(math.pi*math.pi)
values = np.linspace(0,1,f)

fval=[]
for i in range(0,len(values)):
    th = exact(values[i])
    fval.append(th)
plt.plot(values,fval)

fsol = np.vectorize(exact)

#giving boundary conditions

u0 = 0
ub = 1
a_DirNodes = np.array([a, -b])
a_DirVals = np.array([u0, ub])

#--------------------------------------------------------
# use the finite element functions to obtain the solution
#--------------------------------------------------------
conn, f, c = generateMeshNodes(a_Domain, a_Degree, a_Size)


#generate mesh connectivity 

xy = generateMeshConnectivity(conn, a_Degree)


#generate global k matrix

globalstiff = assembleGlobalStiffness(c, xy, a_GaussOrder, a_Degree)
print(globalstiff)

#generate global f matrix

globalload = assembleGlobalLoading(a_Func, c, xy, a_GaussOrder, a_Degree)
print('this is global load')
print(globalload)

#apply dirichlet

[a_Kg, a_Fg] = applyDirichlet(globalstiff, -globalload, a_DirNodes, a_DirVals)

#solve the linear system

[uSol, status] = cg(a_Kg, np.transpose(a_Fg))
print(uSol)



#getting the exact solution to the problem

exactsol = fsol(c)


#compute error 
error = abs(uSol-exactsol)

print('The error for this mesh size is computed as')
print(np.sum(error))



'''plotting the analytical solution'''

#plotting the computed solution
qz = np.array(xy)
#h = plotMesh(, qz)
g = plotSolutions(c, qz, uSol)

