'''
Created on May 17, 2016

@author: carlos
'''

class DataFile(object):
    '''
    classdocs
    '''
    
    def readAscii(self,freq,N_th,N_ph):  
        import numpy as np
        
        import linecache
        
        #use this lines to check if the file is in the folder
        
        '''import os.path
        print(os.listdir(os.getcwd()))
        print(self.filename + '.txt')
        print(os.path.isfile(self.filename + '.txt'))'''
        
        data = []
        for n in range(N_th*N_ph):
            # the "+2" is used to skip the two first lines of the text file
            kk = linecache.getline(self.filename + '.txt', \
                                   3+n+(freq-1)*N_th*N_ph)
            kk = kk.split()
            kk = [np.float(kk[n]) for n in range(len(kk))]
            data = np.append(data, kk, axis=0)
            
        data = np.reshape(data, [N_th*N_ph,8])   
    #    print(data[N_th*N_ph-1,0:4]) 
        self.data = data
        
    def calc_F(self,freq,N_th,N_ph):
        import numpy as np
        
        # reads the data measured at each frequency
        
        self.readAscii(freq,N_th,N_ph)
        
        self.Freq = self.data[0,0]/1e9
        
        self.theta = []
        self.phi = []
        self.E_phi = []
        self.E_the = []
        self.Gain = []
        
        self.theta  = np.array([self.data[n,2]  \
                    for n in range(N_th*N_ph)])
        
        mask = self.data[:,2] < 0
        
        
        
        self.phi  = np.array([self.data[n,1]  \
                    for n in range(N_th*N_ph)])
        
        
        # I change angles from theta -180 to 180 ohi 0 to 180
        # to theta 0 to 180 and phi 0 to 2pi
        # I did this because the integrals where not performed good
        # do not ask me why!
        
        np.copyto(self.phi, np.pi+self.phi, where=mask)
        np.copyto(self.theta, np.abs(self.theta), where=mask)
        
        
        
        self.E_phi = np.array([[self.data[n,3]+1j*self.data[n,4] ] \
                    for n in range(N_th*N_ph)])
        
        self.E_the = np.array([[self.data[n,5]+1j*self.data[n,6] ] \
                    for n in range(N_th*N_ph)])
        
        self.Gain  = np.array([[self.data[n,7] ]  \
                    for n in range(N_th*N_ph)])
        
        
        
        # use this below to test integrals. 
        # Be aware that Gain is converted to linear
        # before integration. UnComment line 22 in Calculations
        
        #self.Gain  = np.array([[self.data[n,7] ] \
        #            for n in range(N_th*N_ph)])

    def write_CST_like(self,work_folder,data,N_antennas,filename):

        text_file = open(work_folder + filename, 'w')
        
        for m in range(N_antennas):
            text_file.write('port' + str(m+1)+ '\n')
            text_file.write('------------------------------------------------'+ '\n')
            
            for n in range(data.shape[0]):
                text_file.write(str(data[n,m,0]) +'\t'+ str(data[n,m,1]) + '\n')
            text_file.write('\n')
        text_file.close()
        
    def writes_RadPat(self,work_folder,filename):
        
        import numpy as np
        
        text_file = open(work_folder + filename, 'w')
        
        text_file.write('Theta [deg.]  Phi   [deg.]  Abs(Grlz)[dB]' \
            + 'Abs(Theta)[dB]  Phase(Theta)[deg.]  Abs(Phi)[dB]  Phase(Phi)[deg.]' +\
            'Ax.Ratio[dB]    \n')
        
        print('not all the info is printed')
        
        text_file.write('------------------------------------------------\n')
            
        for n in range(self.Nphi*self.Ntheta):
            text_file.write(str(np.around(self.theta[n]*180/np.pi, decimals = 2 ))+ '\t' + \
                            str(np.around(self.phi[n]*180/np.pi, decimals = 2 ))+ '\t' + \
                            str(np.around(self.Gain[n], decimals = 4 ))[1:-1] + '\t' + \
                            '0.0000' + '\t' + '0.0000' + '\t' + '0.0000' + '\t' + '0.0000' + '\t' + \
                            '0.0000' + '\n')

        text_file.close()        
    def decimate(self,N_th,N_ph):
        '''
        The data comes from the satimo each 2.5degrees
        I create dummy vectors of same size than the ones coming
        from satimo. Then take only the ones multiple of 15deg
        
        the dummy vectors are created as the come from satimo,
        this is theta -180 to 180 and phi from 0 to 180.
        
        '''
        import numpy as np
        
        
        
        theta = np.array([-180.0+m*360.0/(N_th-1)  \
                    for n in range(N_ph)
                    for m in range(N_th) ])

        phi = np.array([n*180.0/(N_ph)  \
                            for n in range(N_ph)
                            for m in range(N_th) ])
        
        mask1 = np.mod(theta,15) == 0
        
        self.theta = self.theta[mask1]
        self.phi = self.phi[mask1]
        self.E_phi = self.E_phi[mask1]
        self.E_the = self.E_the[mask1]
        self.Gain = self.Gain[mask1]
        
        theta = theta[mask1]
        phi = phi[mask1]
        
        mask2 = np.mod(phi,15) == 0  
        
        self.theta = self.theta[mask2]
        self.phi = self.phi[mask2]
        self.E_phi = self.E_phi[mask2]
        self.E_the = self.E_the[mask2]
        self.Gain = self.Gain[mask2]
        self.Ntheta = 25
        self.Nphi = 12
          
    def __init__(self, filename, Ntheta, Nphi):
        '''
        Constructor
        '''
        self.filename = filename
        self.Ntheta = Ntheta
        self.Nphi = Nphi
    