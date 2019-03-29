'''
Created on May 16, 2016

@author: carlos
'''
import numpy as np

class UnitSphere(object):
    '''
    classdocs
    '''


    def __init__(self, Nth , Nph, theta, phi):
        '''
        Constructor
        '''
#         self.Nt = Nth
#         self.Np = Nph
#         
#         dlt_t= 2*np.pi/Nth
#         dlt_p= np.pi/Nph
#         
#        
#         R = [[np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),\
#               np.cos(theta)] \
#              for theta in np.arange(-np.pi,np.pi,dlt_t) \
#              for phi in np.arange(0,np.pi,dlt_p)]
#     
#         R = np.reshape(R, [np.size(R, 0),3])
        
        R = np.array([[np.sin(theta[n])*np.cos(phi[n]),np.sin(theta[n])*np.sin(phi[n]),\
               np.cos(theta[n])] for n in range(Nth * Nph)])
        
        
        # this are the observation points
        
        self.R = R
        
#         # plot observation points       
#         from mpl_toolkits.mplot3d import Axes3D
#         import matplotlib.pyplot as plt
#           
#         fig = plt.figure()
#         ax = fig.add_subplot(111, projection='3d')
#           
#           
#           
#         ax.set_xlabel('X Label')
#         ax.set_ylabel('Y Label')
#         ax.set_zlabel('Z Label')
#           
#         ax.scatter(R[:,0],R[:,1],R[:,2])
#           
#         plt.show()



        
        from scipy.spatial import ConvexHull

        # compute the convex hull of the points
        cvx = ConvexHull(R)



        dA = []
        for n in range(cvx.nsimplex):
            P0 = cvx.points[cvx.simplices[n,0]]
            P1 = cvx.points[cvx.simplices[n,1]]
            P2 = cvx.points[cvx.simplices[n,2]]
 
            # now i calculate the area of each simplice           
                 
            V1 = np.array(P1 - P0)
            V2 = np.array(P2 - P0)
        
            dA = np.append(dA, np.linalg.norm(np.cross(V1, V2)/2.0))

        self.dA = dA
        self.simplices = cvx.simplices
        self.nsimplex = cvx.nsimplex
