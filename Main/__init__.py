import numpy as np
from telnetlib import theNULL



N_th = 145      
N_ph = 72
Nfreq = 21

Num_ports = 8

import os

# change working directory

work_folder = '/media/carlos/TI31379300A/transfer_VNA/Pollux/Pollux rev 1r03 Set 3/matched  ports/Fields/rawdata/'
processed_folder = '/media/carlos/TI31379300A/transfer_VNA/Pollux/Pollux rev 1r03 Set 3/matched  ports/Fields/processed files/'
os.chdir(work_folder)

# calculate areas of simplices unit Surface for integration


from Satimo import UnitSphere as US
from Satimo import FileManagment as FM

filename = 'port'



# Now I should load the fields at a certain frequency

data1 = FM.DataFile(filename + '1',N_th,N_ph)


# with this I read theta and phi
data1.calc_F(1,N_th,N_ph)


# here I just take the data each 15 degress. 
# I saved data each 2.5deg to make it easy

#data1.decimate(N_th,N_ph)


# with the angles given by SATIMO I calculate the Unit Sphere and all 
#the surface elements dA

sphere1 = US.UnitSphere(data1.Ntheta, data1.Nphi, data1.theta, data1.phi)

print('Test of calculation of area unit sphere')
print('This should be close to  (12.57),', \
      np.round(np.sum(sphere1.dA),2))


import matplotlib.pyplot as plt
from Satimo.Calculations import Calculations 

C1 = Calculations()


eff = np.zeros([Nfreq,Num_ports,2])
rpg = np.zeros([Nfreq,Num_ports,2])

for port in range(Num_ports):
    
    data1 = FM.DataFile(filename + str(port+1),N_th,N_ph)
    
    for f in range(Nfreq):
        
        data1.calc_F(f+1,N_th,N_ph)
        
    #    data1.decimate(N_th,N_ph)
    
    # calculates the Fields and Gain at the center of the simplice
         
        C1.F_at_center(data1, sphere1)
         
    #   C1.calc_ECC(data1, data2, sphere1)
    # calculates efficiency
    
        C1.calc_Eff(data1, sphere1)
    
        
    # calculates peak gain 
        eff[f,port,0] = data1.Freq
        eff[f,port,1] = data1.Eff
    
        rpg[f,port,0] = data1.Freq
        rpg[f,port,1] = np.max(data1.Gain)
    
        
    
        
        if (abs(eff[f,port,0] - 2.45) < 0.0001) :
            
    #        C1.calc_avGphi_Dv(data1,'2G_' + filename1)
    #        C1.calc_avGphi_Dv(data2,'2G_' + filename2)
            data1.writes_RadPat(processed_folder,filename + str(port) + '_Fields@2r4GHz.txt')
    
            
        if (abs(eff[f,port,0] - 5.50) < 0.0001) :
    #       C1.calc_avGphi_Dv(data1,'5G_' + filename1)
    #       C1.calc_avGphi_Dv(data2,'5G_' + filename2)
            
            data1.writes_RadPat(processed_folder,filename + str(port) +'_Fields@5r5GHz.txt')
    
    
                
            
         
        print("%.3f" % eff[f,port,0], "%.3f" %  eff[f,port,1], "%.3f" % rpg[f,port,1])
             
        # save efficiencie0s
        
        
    
    data1.write_CST_like(processed_folder, eff, Num_ports,'eff_2G.txt')  
        
    
    
    data1.write_CST_like(processed_folder, rpg, Num_ports,'RPG_2G.txt')
            