'''
Created on May 19, 2016

@author: carlos
'''

class Calculations(object):
    '''
    classdocs
    '''
    
    def F_at_center(self,data,sphere):
        ''''calculates the field at the center of the simplice'''
        import numpy as np
        
        data.cE_phi = [1.0/3.0*np.sum(data.E_phi[sphere.simplices[k,:]]) \
                         for k in range(sphere.nsimplex)]
        
        data.cE_the = [1.0/3.0*np.sum(data.E_the[sphere.simplices[k,:]]) \
                         for k in range(sphere.nsimplex)]
        
        Gain = np.power(10.0,data.Gain/10.0)
        #Gain = data.Gain
        
        data.cGain = [1.0/3.0*np.sum(Gain[sphere.simplices[k,:]]) \
                         for k in range(sphere.nsimplex)]
        
        

        
    def calc_ECC(self,data1,data2,sphere1):
        import numpy as np 
        
        ''' Please note that the fields are evaluate at the centre of the
         simplices. These Fields are denoted by a c before the components 
         of the field'''
        
        
        # Q12 is the scalar product of E1 and E2 complex conjugate
        # Q1 is the square magnitude of E1 and Q2 the square magnitude of E2 
            
        Q12 = data1.cE_the*np.conjugate(data2.cE_the)+ \
            data1.cE_phi*np.conjugate(data2.cE_phi)
                  
        Q1 = data1.cE_the*np.conjugate(data1.cE_the)+ \
            data1.cE_phi*np.conjugate(data1.cE_phi)
                     
        Q2 = data2.cE_the*np.conjugate(data2.cE_the)+ \
            data2.cE_phi*np.conjugate(data2.cE_phi)           
            
            # Then calculate the integrals
        I12 = np.sum(np.multiply(Q12,sphere1.dA))
        I1 = np.sum(np.multiply(Q1,sphere1.dA))
        I2 = np.sum(np.multiply(Q2,sphere1.dA))
                
        ECC = [np.real(np.abs(I12)/(np.sqrt(I1)*np.sqrt(I2)))]
        
        self.ECC = np.array(ECC)[0]
        
        
                
    def calc_int_rad(self,data):       
        import numpy as np
        
        eta =120.0*np.pi
        data.U = 1.0/(2.0*eta)*(data.E_the*np.conjugate(data.E_the)+ \
                 data.E_phi*np.conjugate(data.E_phi))
        
        
    def calc_Eff(self, data, sphere1):
        import numpy as np
        
        data.Eff = 1.0/(4.0*np.pi)*np.sum(np.multiply(data.cGain,sphere1.dA))
        
    
    def calc_avGphi_Dv(self,data,filename1):    
        '''Calculates some sort of average gain over 
        theta using Diego's approach'''
        import numpy as np
        
        Nth = data.Ntheta
        Nph = data.Nphi
        
        D_th = 2.5*np.pi/180.0
        
        phi_in_deg = np.arange(0,360,2.5)
        av_G_at_phi = np.zeros((np.size(phi_in_deg),))
        # First, select
        k = 0
        for elem in phi_in_deg:
            
            G_at_phi = []
            th_at_phi = []
            shapito = np.shape(data.Gain[ \
                    np.abs(data.phi-elem*np.pi/180.0)<0.001])[0]
                    
            G_at_phi = np.reshape(np.power(10,data.Gain[ \
                    np.abs(data.phi-elem*np.pi/180.0)<0.001]/10.0),shapito)
            th_at_phi = data.theta[ \
                    np.abs(data.phi-elem*np.pi/180.0)<0.001]
            
            av_G_at_phi[k] = 0.5*np.sum(G_at_phi*np.sin(th_at_phi)*D_th)
            k = k + 1
        np.savetxt(filename1+'_Gain_phi.out', [[phi_in_deg[n],av_G_at_phi[n]] \
                    for n in range(np.size(phi_in_deg))], delimiter=',')    
     
     
    def prints_plane_kats(self,data,filename):
#        import matplotlib.pyplot as plt
        import numpy as np
#    extract the data of Gain in plane XZ
        a1 = np.abs(data.phi - np.pi)<0.001
        a2 = np.abs(data.phi - 0)<0.001
       

        theta = np.concatenate((-data.theta[a1],data.theta[a2]),axis = 0)
        Gain= np.concatenate((data.Gain[a1],data.Gain[a2]),axis = 0)
   
        np.savetxt(filename+'_XZ_Gain'+'.out', [[theta[n],Gain[n]] \
            for n in range(np.size(theta))], delimiter=',')
        
#         plt.plot(theta,Gain)
#         plt.show()


#    extract the data of Gain in plane YZ   
        a3 = np.abs(data.phi - 1.5*np.pi)<0.001
        a4 = np.abs(data.phi - 0.5*np.pi)<0.001
        
        theta = np.concatenate((-data.theta[a3],data.theta[a4]),axis = 0)
        Gain= np.concatenate((data.Gain[a3],data.Gain[a4]),axis = 0)        
       
        np.savetxt(filename+'_YZ_Gain'+'.out', [[theta[n],Gain[n]] \
            for n in range(np.size(theta))], delimiter=',')
        
#    extract the data of Gain in plane XY          
        a5 = np.abs(data.theta - 0.5*np.pi)<0.001
        phi = data.phi[a5]
        aux= data.Gain[a5]
        Gain= np.concatenate((Gain[::2],Gain[1::2]),axis = 0)
        
        np.savetxt(filename+'_XY_Gain'+'.out', [[phi[n],Gain[n]] \
                for n in range(np.size(phi))], delimiter=',')      
        
        
        print('232')
    def __init__(self):
        '''
        Constructor
        '''
        