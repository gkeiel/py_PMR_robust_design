# import modules
import math
import numpy as np
import matplotlib.pyplot as plt
from control.matlab import ss, ss2tf 
from scipy.linalg import block_diag
from scipy import signal


def ups_model(r: float, L: float, c: float, R: float, T: float):
    A = [[-(r/L), -1/L],
         [1/c,    -1/(R*c)]]
    B = [[1/L], [0]]
    C = [[0, 1]]
    D = 0
    E = [[0], [-1/c]]

    A = np.eye(2) +T*np.array(A)
    B = T*np.array(B)
    C = np.array(C)
    D = np.array(D)
    E = T*np.array(E)
    return A, B, C, D, E

def pmr_model(f: float, h: int, xi: float, T: float):
    Ar = np.empty((0, 0))
    Br = np.empty((0, 1))
    for i in range(1, h+1):
        wn = 2*math.pi*f*(2*i -1)
        wd = wn*np.sqrt(1 -xi**2)
        AR = [[0,                     1],
              [-math.exp(-2*xi*wn*T), 2*math.cos(wd*T)*math.exp(-xi*wn*T)]]
        BR = [[0], [1]]
        Ar = block_diag(Ar, AR)
        Br = np.vstack((Br, BR))
    return Ar, Br

def agm_model(A, B, C, E, Ar, Br, h: int):
    Aa = np.vstack((np.column_stack((A, np.zeros((2, 2*h)))), np.column_stack((-Br @ C, Ar))))
    Bq = np.vstack((B, Br))
    Ba = np.vstack((B, np.zeros((2*h, 1))))
    Ca = np.column_stack((C, np.zeros((1, 2*h))))
    Da = 0
    Ea = np.vstack((E, np.zeros((2*h, 1))))
    return Aa, Bq, Ba, Ca, Da, Ea

def compute_gain(Aa, Ba, h: int, s: int, T: float):
    n = 2*h +2
    lamb_T = np.linspace(math.exp(-s*T),0.6,n)  # desired eigenvalues
    I = np.eye(n)
    U = np.empty((n, 0))
    W = np.ones(n)
    for i in range(n):
        v = np.linalg.solve(lamb_T[i]*I -Aa, Ba)
        U = np.column_stack((U, v))

    K = np.dot(W, np.linalg.inv(U))
    print("K =", K)
    return K

def get_tf(Aa, Ba, Ca, Da, Ea, Bq, K, h: int, T: float):
    Bk    = Bq*np.vstack((K[1], 0, np.ones((2*h, 1))))
    cl_tf = ss2tf(Aa+Ba*K, Bk, Ca, Da, T)
    id_tf = ss2tf(Aa+Ba*K, Ea, Ca, Da, T)
    return cl_tf, id_tf

def simulate(cl_tf, id_tf, V, R, f, T, Tf=0.2, st=0.04):
    # fundamental sampling time and time vector
    dt = T
    t  = np.linspace(0, Tf-dt, int(Tf/dt))

    # response to sine wave input
    r = V*np.sqrt(2)*np.sin(2*math.pi*f*t)
    r[t <= st] = 0
    t, y1 = signal.dlsim((cl_tf.num[0][0], cl_tf.den[0][0], T), r, t)

    # response to disturbance input
    I_1  = V/R
    I_3  = 0.86*(V/R)
    I_5  = 0.62*(V/R)
    I_7  = 0.35*(V/R)
    I_9  = 0.12*(V/R)
    I_11 = 0.04*(V/R)
    i1   =  I_1*np.sqrt(2)*np.sin(2*math.pi*f*t)
    i3   = -I_3*np.sqrt(2)*np.sin(3*2*math.pi*f*t)
    i5   =  I_5*np.sqrt(2)*np.sin(5*2*math.pi*f*t)
    i7   = -I_7*np.sqrt(2)*np.sin(7*2*math.pi*f*t)
    i9   =  I_9*np.sqrt(2)*np.sin(9*2*math.pi*f*t)
    i11  = -I_11*np.sqrt(2)*np.sin(11*2*math.pi*f*t)
    i1[t <= st] = 0
    i3[t <= st] = 0
    i5[t <= st] = 0
    i7[t <= st] = 0
    i9[t <= st] = 0
    i11[t <= st] = 0
    t, y3 = signal.dlsim((id_tf.num[0][0],id_tf.den[0][0],T), i3, t)
    t, y5 = signal.dlsim((id_tf.num[0][0],id_tf.den[0][0],T), i5, t)
    t, y7 = signal.dlsim((id_tf.num[0][0],id_tf.den[0][0],T), i7, t)
    t, y9 = signal.dlsim((id_tf.num[0][0],id_tf.den[0][0],T), i9, t)
    t, y11 = signal.dlsim((id_tf.num[0][0],id_tf.den[0][0],T), i11, t)

    # output voltage and current for linear load
    v_o = y1
    i_o = i1

    # output voltage and current for non-linear load
    v_o_nl = y1 +y3 +y5 +y7 +y9 +y11
    i_o_nl = i1 +i3 +i5 +i7 +i9 +i11

    return t, r, v_o, i_o, v_o_nl, i_o_nl

def plot_res(t, r, v_o, i_o, v_o_nl, i_o_nl):
    # plot for linear load
    plt.figure(1)
    plt.step(t, r, '--', label='r(kT) [V]')
    plt.step(t, v_o, label='v_o(kT) [V]')
    plt.step(t, i_o, label='i_o(kT) [A]')
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.title('Linear load')
    plt.legend()
    plt.grid(True)

    # plot for non-linear load
    plt.figure(2)
    plt.step(t, r, '--', label='r(kT) [V]')
    plt.step(t, v_o_nl, label='v_o(kT) [V]')
    plt.step(t, i_o_nl, label='i_o(kT) [A]')
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.title('Non-linear load')
    plt.legend()
    plt.grid(True)
    plt.show()
    return 0