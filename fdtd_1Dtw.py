import numpy as np
import h5py
import matplotlib.pyplot as plt

from buildAE4 import buildAE4
from buildAH4 import buildAH4
from buildD2 import buildD2
from buildD343 import buildD343
from buildD4 import buildD4

np.set_printoptions(precision=12)
np.set_printoptions(suppress=True)
np.set_printoptions(linewidth=150)
tol = 1e-12


def sinewave(x, t):
    return np.sin(3 * np.pi * (x+t))


dx=1/40
dt=1/1600
Tmax=4



def fdtd1Dtw(method, dx=1/40, dt=1/1600, Tmax=4):

    nmax=int(np.ceil(Tmax/dt))
    I = int(np.ceil(1/dx))

    Hy = np.zeros(I)
    Ez = np.zeros(I+1)

    errEz = np.zeros(nmax)
    L2errEz = np.zeros(nmax)

    if method == 'cg3434':
        # CG:4/343
	Deh, Q = buildD4(I)
	Dhe, Q = buildD343(I-1);
    elif method == 'explicit24':
        # exact(2,4)
        Deh = buildAE4(I)
        Dhe = buildAH4(I)
    elif method == 'cg4':
        # CG:4/4
        Deh, Q = buildD4(I);
        Dhe, Q  = buildD4(I-1)
    elif method == 'cg343':
        # CG:343/343
        Deh, Q = buildD343(I)
        Dhe, Q = buildD343(I-1)
    elif method == 'yee222':
        # standard FDTD
        dt=dx
        Deh = buildD2(I)
        Dhe = buildD2(I-1)


    n = 0
    # Ez @ t = 0
    EzIC = sinewave(np.arange(1, I)*dx, n*dt)
    

    n = .5
    # Hy @ t = (1/2)dt
    HyIC = sinewave((np.arange(0, I)+.5)*dx, n*dt)
    
    n = 1
    Ez[1:I] = EzIC + (dt/dx) * Dhe * HyIC
    
    # BCs
    Ez[0] = sinewave(0, n*dt);
    Ez[I] = sinewave(1, n*dt);
    
    errEz = np.abs(Ez-sinewave(np.arange(I+1)*dx,n*dt))
    L2errEz[0] = np.sqrt(dx)*np.linalg.norm(errEz,2)

    # Hy @ t=(3/2)dt
    Hy = HyIC + (dt/dx) * Deh * Ez

    for n in np.arange(2,nmax+1):
        # Ez @ t = n*dt
        Ez[1:I] = Ez[1:I] + (dt/dx) * (Dhe * Hy)
        # BCs
        Ez[0] = sinewave(0,n*dt);
        Ez[I] = sinewave(1,n*dt);
        #Hy @ t = (n+1/2)dt
        Hy = Hy + (dt/dx) * Deh * Ez;
        #Error
        errEz = np.abs(Ez-sinewave(np.arange(I+1)*dx,n*dt))
        L2errEz[n-1] = np.sqrt(dx)*np.linalg.norm(errEz,2)


    # COMPARISON
    #data = h5py.File('L2errEz.arr', 'r')
    #octL2errEz = data['L2errEz']['value'].value.ravel()

    #data = h5py.File('Ez.arr', 'r')
    #octEz =  data['Ez']['value'].value.ravel()

    #data = h5py.File('Hy.arr', 'r')
    #octHy =  data['Hy']['value'].value.ravel()

    #data = h5py.File('fDeh.mat', 'r')
    #octDeh =  data['fDeh']['value'].value.T
    #data.close()

    #data = h5py.File('fDhe.mat', 'r')
    #octDhe =  data['fDhe']['value'].value.T
    #data.close()

    # np.all(abs(Dhe - octDhe) < tol)
    # np.all(abs(Deh - octDeh) < tol)
    # np.all(abs(Ez - octEz) < tol)
    # np.all(abs(L2errEz - octL2errEz) < tol)


    # # PLOTS
    # plt.plot(Ez)
    # plt.plot(errEz)
    # plt.plot(L2errEz)
    # plt.plot(abs(L2errEz-octL2errEz))
    # plt.plot(QQ[1:])
    # plt.ylim((0,0.05))
    # plt.show()
