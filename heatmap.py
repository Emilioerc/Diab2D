import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import os


''' ADC Surfaces '''

GS     = []
ES1    = []
ES2    = []
ES3    = []
ES4    = []
ES5    = []
r      = []
theta  = []

data0  = np.loadtxt('theta2.txt')
data1  = np.loadtxt('r2.txt')
data2  = np.loadtxt('GS_pes.txt')
data3  = np.loadtxt('ES1_pes.txt')
data4  = np.loadtxt('ES2_pes.txt')
data5  = np.loadtxt('ES3_pes.txt')
data6  = np.loadtxt('ES4_pes.txt')
data7  = np.loadtxt('ES5_pes.txt')


theta.append(data0)
r.append(data1)
GS.append(data2)
ES1.append(data3)
ES2.append(data4)
ES3.append(data5)
ES4.append(data6)
ES5.append(data7)

#energy = np.array(energy)
#energy = np.abs(energy)
#print('ES1 max:', np.sort(energy))

tt, rr = np.meshgrid(theta, r)
GS  = np.array(GS)
ES1 = np.array(ES1)
ES2 = np.array(ES2)
ES3 = np.array(ES3)
ES4 = np.array(ES4)
ES5 = np.array(ES5)
GS  = GS.reshape(38,10) # first dim runs over r, second over theta
ES1 = ES1.reshape(38,10)
ES2 = ES2.reshape(38,10)
ES3 = ES3.reshape(38,10)
ES4 = ES4.reshape(38,10)
ES5 = ES5.reshape(38,10)

print('rr', rr)
print('tt', tt)
print('GS', GS)


''' CASSCF Surfaces, Faraday Discuss., 2004, 127, 283-293 '''


def V11(r, theta, p1, p2, p3, d1, d2, beta1, alpha1, De, a1, r1):
    f1 = 0.5 * (1 + np.tanh((r - d2) / beta1))
    k1 = (p1 + p2 * r) * (1 - f1) + p3 * np.exp(-(r - d1) / alpha1) * f1
    V11 = De * (1 - np.exp(-a1 * (r - r1)))**2 + 0.5 * k1 * theta**2
    return V11

def V22(r, theta, E02, De2in, r2, a2, De2, a3, r3, alpha2, q1, q2, q3, q4, l22):
    Vin = E02 + De2in*(1-np.exp(-a2*(r-r2)))**2
    Vout = De2 + a3*np.exp(-(r-r3)/alpha2)
    k2 = 0.5*(q1+q2*r) - 0.5*np.sqrt((q3 + q2*r)**2 +4*q4**2)
    k2 = np.where(r>2.55, 0, k2)
    V22 = 0.5*(Vin + Vout) - 0.5*np.sqrt((Vin - Vout)**2 + 4*l22**2) + 0.5*k2*theta**2
    return V22

r = np.array(r)
theta = np.array(theta)


#Parameters V11
De = 5.117 #eV
r1 = 1.959 #au
a1 = 1.196 #au
p1 = 5.147 #eV
p2 = -1.344 #eV au^-1
p3 = 0.884 #eV
alpha1 = 0.775 #au
beta1  = 0.00015 #au
d1     = 3.100 #au
d2     = 2.696 #au

#Parameters V22
E02 = 5.584 
De2in = 8.070
r2 = 1.922
a2 = 0.882
De2 = 4.092
a3 = 0.091
r3 = 5.203
alpha2 = 0.774
q1 = 3.818
q2 = -1.219
q3 = 2.444
q4 = 0.226
l22 = 1.669


r = np.linspace(1.7,10,100)
theta = np.linspace(-2, 2, 100)
r, theta = np.meshgrid(r, theta)


V_GS   = V11(r, theta, p1, p2, p3, d1, d2, beta1, alpha1, De, a1, r1)
V_ES   = V22(r, theta, E02, De2in, r2, a2, De2, a3, r3, alpha2, q1, q2, q3, q4, l22) 

print('r:',r)
print('theta:', theta)
print('V11:', V11)


fig, ax = plt.subplots(subplot_kw={"projection":"3d"})


surf = ax.plot_surface(r,theta, V_GS, cmap=cm.coolwarm, linewidth=0,antialiased=False)
surf2= ax.plot_surface(r,theta, V_ES, cmap=cm.coolwarm, linewidth=0,antialiased=False)
plt.show()





''' FIGURES '''


fig, ax = plt.subplots(subplot_kw={"projection":"3d"})


surf = ax.plot_surface(tt,rr,GS, cmap=cm.coolwarm, linewidth=0,antialiased=False)
surf2= ax.plot_surface(tt,rr,ES1,cmap=cm.coolwarm, linewidth=0,antialiased=False)
surf3= ax.plot_surface(tt,rr,ES2,cmap=cm.coolwarm, linewidth=0,antialiased=False)
#surf4= ax.plot_surface(tt,rr,ES3,cmap=cm.coolwarm, linewidth=0,antialiased=False)
#surf5= ax.plot_surface(tt,rr,ES4,cmap=cm.coolwarm, linewidth=0,antialiased=False)
#surf6= ax.plot_surface(tt,rr,ES5,cmap=cm.coolwarm, linewidth=0,antialiased=False)

plt.show()

exit()

plt.figure(figsize=(8,8))
plt.subplot(2,1,1)
plt.contourf(rr, tt, GS, cmap='RdBu',levels=50)
plt.xlabel('r')
plt.ylabel('theta')
plt.colorbar()
plt.title('Ground State')

plt.subplot(2,1,2)
plt.contourf(rr[:,:], tt[:,:], ES1[:,:], cmap='RdBu',levels=50)
plt.xlabel('r')
plt.ylabel('theta')
plt.colorbar()
plt.title('1st Excited State')
plt.show()


exit()


# Create the heatmap
plt.figure(figsize=(10, 6))
plt.subplot(3,4,1)
plt.contourf(tt[5,:,:], dd[5,:,:], energy[5,:,:], cmap='RdBu',vmin=0,vmax=0.7 ,levels=50)
plt.xlabel('theta')
plt.ylabel('rd')
plt.title('rv=2.1')

plt.subplot(3,4,2)
plt.contourf(tt[6,:,:], dd[6,:,:], energy[6,:,:], cmap='RdBu',vmin=0,vmax=0.7 ,levels=50)
plt.xlabel('theta')
plt.ylabel('rd')
plt.title('rv=2.2')

plt.subplot(3,4,3)
plt.contourf(tt[7,:,:], dd[7,:,:], energy[7,:,:], cmap='RdBu',vmin=0,vmax=0.7 ,levels=50)
plt.xlabel('theta')
plt.ylabel('rd')
plt.title('rv=2.3')

plt.subplot(3,4,4)
plt.contourf(tt[8,:,:], dd[8,:,:], energy[8,:,:], cmap='RdBu',vmin=0,vmax=0.7 ,levels=50)
plt.xlabel('theta')
plt.ylabel('rd')
#plt.colorbar(label='Energy (au)')
plt.title('rv=2.4')


plt.subplot(3,4,5)
plt.contourf(tt[:,2,:], vv[:,2,:], energy[:,2,:], cmap='RdBu',vmin=0,vmax=0.7 ,levels=50)
plt.xlabel('theta')
plt.ylabel('rv')
plt.title('rd=4.24')


plt.subplot(3,4,6)
plt.contourf(tt[:,3,:], vv[:,3,:], energy[:,3,:], cmap='RdBu',vmin=0,vmax=0.7 ,levels=50)
plt.xlabel('theta')
plt.ylabel('rv')
plt.title('rd=4.46')

plt.subplot(3,4,7)
plt.contourf(tt[:,4,:], vv[:,4,:], energy[:,4,:], cmap='RdBu',vmin=0,vmax=0.7 ,levels=50)
plt.xlabel('theta')
plt.ylabel('rv')
plt.title('rd=4.68')

plt.subplot(3,4,8)
plt.contourf(tt[:,5,:], vv[:,5,:], energy[:,5,:], cmap='RdBu',vmin=0,vmax=0.7 ,levels=50)
plt.xlabel('theta')
plt.ylabel('rv')
plt.title('rd=4.9')
#plt.colorbar(label='Energy (au)')


plt.subplot(3,4,9)
plt.contourf(vv[:,:,4], dd[:,:,4], energy[:,:,4], cmap='RdBu',vmin=0,vmax=0.7,  levels=50)
plt.xlabel('rv')
plt.ylabel('rd')
plt.title('theta=1.57')

plt.subplot(3,4,10)
plt.contourf(vv[:,:,5], dd[:,:,5], energy[:,:,5], cmap='RdBu',vmin=0,vmax=0.7,  levels=50)
plt.xlabel('rv')
plt.ylabel('rd')
plt.title('theta=1.25')

plt.subplot(3,4,11)
plt.contourf(vv[:,:,6], dd[:,:,6], energy[:,:,6], cmap='RdBu',vmin=0,vmax=0.7,  levels=50)
plt.xlabel('rv')
plt.ylabel('rd')
plt.title('theta=0.94')

plt.subplot(3,4,12)
plt.contourf(vv[:,:,7], dd[:,:,7], energy[:,:,7], cmap='RdBu',vmin=0,vmax=0.7,  levels=50)
plt.xlabel('rv')
plt.ylabel('rd')
#plt.colorbar(label='Energy (au)')
plt.title('theta=0.67')


#plt.show()



# Convert the absorption data to a numpy array
#absorption_data = np.array(absorption_data)
# Assuming that energy values are the same for all files
#energy_values = data[:, 0]



n = 10
#colors = plt.cm.RdBu(np.linspace(0.1,n))



plt.figure(figsize=(10,6))
plt.subplot(1,5,1)
plt.plot(np.reshape(rd, (10,)), energy[4,:,0], label='theta=max')# color=colors[0])
plt.plot(np.reshape(rd, (10,)), energy[4,:,1]) #color=colors[1])
plt.plot(np.reshape(rd, (10,)), energy[4,:,2]) #color=colors[2])
plt.plot(np.reshape(rd, (10,)), energy[4,:,3]) #color=colors[3])
plt.plot(np.reshape(rd, (10,)), energy[4,:,4]) #color=colors[4])
plt.plot(np.reshape(rd, (10,)), energy[4,:,5]) #color=colors[5])
plt.plot(np.reshape(rd, (10,)), energy[4,:,6]) #color=colors[6])
plt.plot(np.reshape(rd, (10,)), energy[4,:,7]) #color=colors[7])
plt.plot(np.reshape(rd, (10,)), energy[4,:,8]) #color=colors[8])
plt.plot(np.reshape(rd, (10,)), energy[4,:,9], label='theta=min') # color=colors[9])
plt.title('rv=2')
plt.legend()
plt.xlabel('rd')

plt.subplot(1,5,2)
plt.plot(np.reshape(rd, (10,)), energy[5,:,0])
plt.plot(np.reshape(rd, (10,)), energy[5,:,1])
plt.plot(np.reshape(rd, (10,)), energy[5,:,2])
plt.plot(np.reshape(rd, (10,)), energy[5,:,3])
plt.plot(np.reshape(rd, (10,)), energy[5,:,4])
plt.plot(np.reshape(rd, (10,)), energy[5,:,5])
plt.plot(np.reshape(rd, (10,)), energy[5,:,6])
plt.plot(np.reshape(rd, (10,)), energy[5,:,7])
plt.plot(np.reshape(rd, (10,)), energy[5,:,8])
plt.plot(np.reshape(rd, (10,)), energy[5,:,9])
plt.title('rv=2.1')
plt.xlabel('rd')

plt.subplot(1,5,3)
plt.plot(np.reshape(rd, (10,)), energy[6,:,0])
plt.plot(np.reshape(rd, (10,)), energy[6,:,1])
plt.plot(np.reshape(rd, (10,)), energy[6,:,2])
plt.plot(np.reshape(rd, (10,)), energy[6,:,3])
plt.plot(np.reshape(rd, (10,)), energy[6,:,4])
plt.plot(np.reshape(rd, (10,)), energy[6,:,5])
plt.plot(np.reshape(rd, (10,)), energy[6,:,6])
plt.plot(np.reshape(rd, (10,)), energy[6,:,7])
plt.plot(np.reshape(rd, (10,)), energy[6,:,8])
plt.plot(np.reshape(rd, (10,)), energy[6,:,9])
plt.title('rv=2.2')
plt.xlabel('rd')

plt.subplot(1,5,4)
plt.plot(np.reshape(rd, (10,)), energy[7,:,0])
plt.plot(np.reshape(rd, (10,)), energy[7,:,1])
plt.plot(np.reshape(rd, (10,)), energy[7,:,2])
plt.plot(np.reshape(rd, (10,)), energy[7,:,3])
plt.plot(np.reshape(rd, (10,)), energy[7,:,4])
plt.plot(np.reshape(rd, (10,)), energy[7,:,5])
plt.plot(np.reshape(rd, (10,)), energy[7,:,6])
plt.plot(np.reshape(rd, (10,)), energy[7,:,7])
plt.plot(np.reshape(rd, (10,)), energy[7,:,8])
plt.plot(np.reshape(rd, (10,)), energy[7,:,9])
plt.title('rv=2.3')

plt.subplot(1,5,5)
plt.plot(np.reshape(rd, (10,)), energy[8,:,0])
plt.plot(np.reshape(rd, (10,)), energy[8,:,1])
plt.plot(np.reshape(rd, (10,)), energy[8,:,2])
plt.plot(np.reshape(rd, (10,)), energy[8,:,3])
plt.plot(np.reshape(rd, (10,)), energy[8,:,4])
plt.plot(np.reshape(rd, (10,)), energy[8,:,5])
plt.plot(np.reshape(rd, (10,)), energy[8,:,6])
plt.plot(np.reshape(rd, (10,)), energy[8,:,7])
plt.plot(np.reshape(rd, (10,)), energy[8,:,8])
plt.plot(np.reshape(rd, (10,)), energy[8,:,9])
plt.title('rv=2.4')




#plt.show()



# Create the heatmap
plt.figure(figsize=(20, 6))
plt.subplot(2,5,1)
plt.contourf(vv[:,:,0], dd[:,:,0], energy[:,:,0], cmap='RdBu',vmin=0,vmax=0.7, levels=50)
plt.xlabel('rv')
plt.ylabel('rd')
plt.title('theta=2.83')

plt.subplot(2,5,2)
plt.contourf(vv[:,:,1], dd[:,:,1], energy[:,:,1], cmap='RdBu',vmin=0,vmax=0.7, levels=50)
plt.xlabel('rv')
plt.ylabel('rd')
plt.title('theta=2.51')

plt.subplot(2,5,3)
plt.contourf(vv[:,:,2], dd[:,:,2], energy[:,:,2], cmap='RdBu',vmin=0,vmax=0.7,  levels=50)
plt.xlabel('rv')
plt.ylabel('rd')
plt.title('theta=2.2')

plt.subplot(2,5,4)
plt.contourf(vv[:,:,3], dd[:,:,3], energy[:,:,3], cmap='RdBu',vmin=0,vmax=0.7,  levels=50)
plt.xlabel('rv')
plt.ylabel('rd')
plt.title('theta=1.88')

plt.subplot(2,5,5)
plt.contourf(vv[:,:,4], dd[:,:,4], energy[:,:,4], cmap='RdBu',vmin=0,vmax=0.7,  levels=50)
plt.xlabel('rv')
plt.ylabel('rd')
plt.colorbar(label='Energy (au)')
plt.title('theta=1.57')

plt.subplot(2,5,6)
plt.contourf(vv[:,:,5], dd[:,:,5], energy[:,:,5], cmap='RdBu',vmin=0,vmax=0.7,  levels=50)
plt.xlabel('rv')
plt.ylabel('rd')
plt.title('theta=1.25')

plt.subplot(2,5,7)
plt.contourf(vv[:,:,6], dd[:,:,6], energy[:,:,6], cmap='RdBu',vmin=0,vmax=0.7,  levels=50)
plt.xlabel('rv')
plt.ylabel('rd')
plt.title('theta=0.94')

plt.subplot(2,5,8)
plt.contourf(vv[:,:,7], dd[:,:,7], energy[:,:,7], cmap='RdBu',vmin=0,vmax=0.7,  levels=50)
plt.xlabel('rv')
plt.ylabel('rd')
plt.title('theta=0.67')


plt.subplot(2,5,9)
plt.contourf(vv[:,:,8], dd[:,:,8], energy[:,:,8], cmap='RdBu',vmin=0,vmax=0.7,  levels=50)
plt.xlabel('rv')
plt.ylabel('rd')
plt.title('theta=0.31')

plt.subplot(2,5,10)
plt.contourf(vv[:,:,9], dd[:,:,9], energy[:,:,9], cmap='RdBu',vmin=0,vmax=0.7,  levels=50)
plt.xlabel('rv')
plt.ylabel('rd')
plt.title('theta=0.0')




plt.show()

