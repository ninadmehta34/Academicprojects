import sys

try:
    import numpy as np
    import matplotlib.pyplot as plt
except ImportError:
    sys.exit('Scientific Python Stack: NumPy+SciPy+MatPlotLib Not Installed!')

try:
    from FEM_1D_Functions import *
except ImportError:
    sys.exit('Functions From FEM_1D_Functions Not Successfully Loaded!')
             
             
a_Degree = int(input("Enter the polynomial order: "))
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
             


# In[ ]:




