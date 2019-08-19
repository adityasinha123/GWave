import numpy as np
import matplotlib.pyplot as plt

#Constants (in SI)
G = 6.6*(10**(-11))#SI
c = 3*(10**8)#SI

#-----------------------Functions-------------#
def binarymass(m1,m2):
    mu =(2*(10**30)) * m1*m2/(m1+m2)#In SI
    M_tot = (2*(10**30)) * m1+m2#in SI
    M_ =(mu)**(3/5) * (M_tot)**(2/5)#((m1*m2)/(m1+m2))**3/5 * (m1 + m2)**2/5 #Chirp mass in SI
    K = G*M_/(c**2)#[K]=L
    return mu,M_tot,M_,K

def BeamPatternFunc(theta,phi,psi):
    #Location of source (r,theta,phi)
    #psi is the Euler angle that translates from detector frame to radiation frame
    F_plus = (-1/2)*(1 + np.cos(theta)**2)*np.cos(2*phi)*np.cos(2*psi) - np.cos(theta)*np.cos(2*phi)*np.sin(2*psi)
    F_x = (1/2)*(1 + np.cos(theta)**2)*np.cos(2*phi)*np.sin(2*psi) - np.cos(theta)*np.sin(2*phi)*np.cos(2*psi)
    return F_plus,F_x

def strain_signal(F_1,F_2,h_1,h_2):
    return F_1*h_1 + F_2*h_2

#----------------------------------------------#
#Mass
mu,M_tot,M_,K = binarymass(1,1)#Mass parameters

#Distance and angle
r = (3.086*10**(22)) * (0.0010)#in m; () in Mega Parsec
iota = np.pi/2#Angle in radian between orbital plane(through y-axis) and z axis
A = (4 * (K)**(5/3))/r#Include the cosmo redshift effect

#time
T,step = np.linspace(0,60,10000, retstep=True)#duration and timse interval in seconds
t = -1*(10**7)+T #Looking at an interval of T at t[0]. 

tc = 0#max(t)+step#tc = 0 coalescence at t=0. Then looking at time interval T, t[0] time before tc(coalescence)

#Frequency and phase
theta =  (mu/M_tot)*((c**3)/(5*G*M_tot))*(tc - t)
phi_c = 0#np.pi/2 #(phase at tc)

f = ((c**3)/(8*np.pi*G*M_tot))*(theta**(-3/8))
phi = (-2*c**3 / (5*G*M_tot))*(theta**(-3/8))*(tc-t) + phi_c


#2PN phi
Phi2 = phi_c - (1/(mu/M_tot))*(theta**(5/8) + (3715/8064 + (mu/M_tot)*55/96)*(theta**(3/8)) - 3*np.pi/4 *theta**(1/4)*(9275495/14450688 + 284875/258048*(mu/M_tot) + (mu/M_tot)*(mu/M_tot)*1855/2048)*theta**(1/8))
f2 = ((c**3)/(8*np.pi*G*M_tot))*(theta**(-3/8) + (743/2688 + 11/32*(mu/M_tot))*theta**(-5/8) - 3*np.pi/10*theta**(-3/4) + (1855099/14450688 + (mu/M_tot)*56975/258048 + 371/2048*((mu/M_tot)**2))*theta**(-7/8))


#Polarizatio(mu/M_tot)s
h_plus = A*((np.pi*f/c)**(2/3))*np.cos(phi)*((1+np.cos(iota)**2)/2)

h_x = A*((np.pi*f/c)**(2/3))*np.sin(phi)*np.cos(iota)

#Beam Pattern
F_plus,F_x = BeamPatternFunc(0,0,0)#At zenith 

h = strain_signal(F_plus,F_x, h_plus,h_x)


#Plots
plt.plot(t,h)
#plt.plot(t,h_x,'r')
plt.title("Gravitational wave Polarizations")
plt.show()


"""
choosing such that the collision happens at t=0 and we are ooking at few weeks before that
3.154*10**7 sec = 1 year,
604800 sec = 1 week
86400 sec = 1 day
"""
