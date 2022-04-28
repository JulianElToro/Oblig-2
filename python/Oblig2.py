import numpy as np
from typing import List
import matplotlib.pyplot as plt

#%% A)

# T/J = 0.25

#The data gets extracted from a .txt file

with  open("Correlation_length_T025.txt", "r") as  infile:

    lines = infile.readlines()

    r: List[float] = []
    C: List[float] = []

    for  line in  lines:
        vals = line.split()
        r.append(float(vals[0]))
        C.append(float(vals[1]))
        

# Mean values for the different measurements bins   
N = 16

Nbins = 10

Cmean = np.zeros(N)

Ci = np.zeros(Nbins)

for i in range(0,N):
    for j in range(0,Nbins):
        
        Ci[j] = C[i+j*N] 
        #print(Ci[j], i+j*N)
            
    Cmean[i] = sum(Ci)/Nbins
    
#Theoretical values

J = 1

T = 0.25 #temperature in units of J/k_b

L1 = np.exp(J/T) - 1 #first eigenvalue

L2 = np.exp(J/T) + 2 #second eigenvalue

r1 = np.linspace(0 , N-1 , N)

Ct_Nf = np.zeros(N)

for i in range(0,N): #Correlation length finite N
    
    Ct_Nf[i] =  (L1**r1[i] * L2**(N-r1[i]) + L1**(N-r1[i]) * L2**r1[i] + L1**N) / (2*L1**N + L2**N)

Ct_Ninf =  ( L1 / L2 ) ** r1 #Correlation length infinite N

#The data and theoretical values get ploted

plt.figure()

plt.plot(r1, Cmean , 'b', label='Simulation', linewidth=5, alpha=.5)
plt.plot(r1 , Ct_Nf, 'r', label='Theoretical for finite $N$')
plt.plot(r1, Ct_Ninf, 'y', label='Theoretical for $N=\infty$')
plt.title('Correlation length (real part) as a function of $r$ for $T/J = 0.25$')
plt.xlabel(r'r')
plt.ylabel(r'C(r)')
plt.legend()
plt.grid(True)

plt.savefig('Oblig_Fig1.pdf')

# T/J = 0.5

#The data gets extracted from a .txt file

with  open("Correlation_length_T05.txt", "r") as  infile:

    lines = infile.readlines()

    r: List[float] = []
    C: List[float] = []

    for  line in  lines:
        vals = line.split()
        r.append(float(vals[0]))
        C.append(float(vals[1]))
        

# Mean values for the different measurements bins  
    
N = 16

Nbins = 10

Cmean = np.zeros(N)

Ci = np.zeros(Nbins)

for i in range(0,N):
    for j in range(0,Nbins):
        
        Ci[j] = C[i+j*N] 
        #print(Ci[j], i+j*N)
            
    Cmean[i] = sum(Ci)/Nbins

# Theoretical values

J = 1

T = 0.5 #temperature in units of J/k_b

L1 = np.exp(J/T) - 1 #first eigenvalue

L2 = np.exp(J/T) + 2 #second eigenvalue

r1 = np.linspace(0 , N-1 , N)

Ct_Nf = np.zeros(N)

for i in range(0,N): # Correlation length finite N
    
    Ct_Nf[i] =  (L1**r1[i] * L2**(N-r1[i]) + L1**(N-r1[i]) * L2**r1[i] + L1**N) / (2*L1**N + L2**N)

Ct_Ninf =  ( L1 / L2 ) ** r1 # Correlation length infinite N

#The data and theretical values get ploted

plt.figure()

plt.plot(r1, Cmean , 'b', label='Simulation' , linewidth=5, alpha=.5)
plt.plot(r1 , Ct_Nf, 'r', label='Theoretical for finite $N$')
plt.plot(r1, Ct_Ninf, 'y', label='Theoretical for $N=\infty$')
plt.title('Correlation length (real part) as a function of $r$ for $T/J = 0.5$')
plt.xlabel(r'r')
plt.ylabel(r'C(r)')
plt.legend()
plt.grid(True)

plt.savefig('Oblig_Fig2.pdf')


# %% C)

#The data gets extracted from a .txt file
       
with  open("mvsTJ_L16_b.txt", "r") as  infile:

    lines = infile.readlines()

    T_16_b: List[float] = [] #temperature in units of J/k_b
    mRe_16_b: List[float] = [] #average magnetization
    m_abs_16_b: List[float] = [] #average absolute magnetization
    m2_16_b: List[float] = [] #average absolute magnetization squared
    m4_16_b: List[float] = [] #average absolute magnetization to the 4th

    for  line in  lines:
        vals = line.split()
        T_16_b.append(float(vals[0]))
        mRe_16_b.append(float(vals[1]))
        m_abs_16_b.append(float(vals[2]))
        m2_16_b.append(float(vals[3]))
        m4_16_b.append(float(vals[4]))
        
   
#The data gets ploted

plt.figure()

plt.plot(T_16_b, mRe_16_b , 'b', label='Real part of m')
plt.title(r'Real part of the magnetization $m$ as a function of $T/J$')
plt.xlabel(r'$T/J$')
plt.ylabel(r'$\Re(m)$')
plt.grid(True)

plt.savefig('Oblig_Fig3_b.pdf')

#%% D)

plt.figure()

plt.plot(T_16_b, m2_16_b, 'b', label='average absolute magnetization squared')
plt.title(r'Average absolute magnetization squared $\langle |m|^2 \rangle$ as a function of $T/J$')
plt.xlabel(r'$T/J$')
plt.ylabel(r'$\langle |m|^2 \rangle$')
plt.grid(True)

plt.savefig('Oblig_Fig4_b.pdf')
      

    
# %% F)

#The data gets extracted from different .txt files

with  open("mvsTJ_L8.txt", "r") as  infile:

    lines = infile.readlines()

    T_L8: List[float] = [] #temperature in units of J/k_b
    mRe_L8: List[float] = [] #average magnetization
    m_abs_L8: List[float] = [] #average absolute magnetization
    m2_L8: List[float] = [] #average absolute magnetization squared
    m4_L8: List[float] = [] #average absolute magnetization to the 4th

    for  line in  lines:
        vals = line.split()
        T_L8.append(float(vals[0]))
        mRe_L8.append(float(vals[1]))
        m_abs_L8.append(float(vals[2]))
        m2_L8.append(float(vals[3]))
        m4_L8.append(float(vals[4]))

T_L8 = np.asarray(T_L8)
m2_L8 = np.asarray(m2_L8)
m4_L8 = np.asarray(m4_L8)

with  open("mvsTJ_L16.txt", "r") as  infile:

    lines = infile.readlines()

    T_L16: List[float] = [] #temperature in units of J/k_b
    mRe_L16: List[float] = [] #average magnetization
    m_abs_L16: List[float] = [] #average absolute magnetization
    m2_L16: List[float] = [] #average absolute magnetization squared
    m4_L16: List[float] = [] #average absolute magnetization to the 4th

    for  line in  lines:
        vals = line.split()
        T_L16.append(float(vals[0]))
        mRe_L16.append(float(vals[1]))
        m_abs_L16.append(float(vals[2]))
        m2_L16.append(float(vals[3]))
        m4_L16.append(float(vals[4]))        
        
T_L16 = np.asarray(T_L16)
m2_L16 = np.asarray(m2_L16)
m4_L16 = np.asarray(m4_L16)


with  open("mvsTJ_L32.txt", "r") as  infile:

    lines = infile.readlines()

    T_L32: List[float] = [] #temperature in units of J/k_b
    mRe_L32: List[float] = [] #average magnetization
    m_abs_L32: List[float] = [] #average absolute magnetization
    m2_L32: List[float] = [] #average absolute magnetization squared
    m4_L32: List[float] = [] #average absolute magnetization to the 4th

    for  line in  lines:
        vals = line.split()
        T_L32.append(float(vals[0]))
        mRe_L32.append(float(vals[1]))
        m_abs_L32.append(float(vals[2]))
        m2_L32.append(float(vals[3]))
        m4_L32.append(float(vals[4]))

T_L32 = np.asarray(T_L32)
m2_L23 = np.asarray(m2_L32)
m4_L32 = np.asarray(m4_L32)

#Gamma

len(T_L8)

G_L8 = np.zeros(len(T_L8))
G_L16 = np.zeros(len(T_L16))
G_L32 = np.zeros(len(T_L32))

for i in range(0 , len(T_L8)):
    G_L8[i] = m4_L8[i] /( m2_L8[i]**2)
    G_L16[i] =   m4_L16[i] /(m2_L16[i]**2)
    
for i in range(0 , len(T_L32)):
    G_L32[i] = m4_L32[i] /(m2_L32[i]**2)
  
#Plot

plt.figure()

plt.plot(T_L8, G_L8 , 'b', label=r'$\Gamma$ for $L=8$')
plt.plot(T_L16, G_L16 , 'r', label=r'$\Gamma$ for $L=16$')
plt.plot(T_L32, G_L32, 'g', label=r'$\Gamma$ for $L=32$')

plt.title(r'$\Gamma $ as a function of $T/J$ for different $L$')
plt.xlabel(r'$T/J$')
plt.ylabel(r'$\Gamma$')
plt.legend()
plt.grid(True)

plt.savefig('Oblig_Fig5.pdf')
        
        