# import modules
import math
import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt
from control.matlab import ss, ss2tf 
from scipy.linalg import block_diag
from scipy import signal


def ups_model(r: float, L: float, c: float, R: float):
    A1 = [[-(r/L), -1/L],
          [1/c,     0]]
    A2 = [[-(r/L), -1/L],
          [1/c,    -1/(R*c)]]
    B = [[1/L], [0]]
    C = [[0, 1]]
    D = 0
    E = [[0], [-1/c]]
    return np.array(A1), np.array(A2), np.array(B), np.array(C), np.array(D), np.array(E)

def pmr_model(f: float, h: int, xi: float):
    Ar = np.empty((0, 0))
    Br = np.empty((0, 1))
    for i in range(1, h+1):
        wn = 2*math.pi*f*(2*i -1)
        AR = [[0,    wn],
              [-wn, -2*xi*wn]]
        BR = [[0], [1]]
        Ar = block_diag(Ar, AR)
        Br = np.vstack((Br, BR))
    return Ar, Br

def agm_model(A1, A2, B, C, E, Ar, Br, h: int):
    Aa1 = np.vstack((np.column_stack((A1, np.zeros((2, 2*h)))), np.column_stack((-Br @ C, Ar))))
    Aa2 = np.vstack((np.column_stack((A2, np.zeros((2, 2*h)))), np.column_stack((-Br @ C, Ar))))
    Bq = np.vstack((B, Br))
    Ba = np.vstack((B, np.zeros((2*h, 1))))
    Ca = np.column_stack((C, np.zeros((1, 2*h))))
    Da = 0
    Ea = np.vstack((E, np.zeros((2*h, 1))))
    return Aa1, Aa2, Bq, Ba, Ca, Da, Ea

def compute_gain(Aa1, Aa2, Ba, h: int, s: int):
    n = 2*h+2
    p  = 1
    Q  = cp.Variable((n,n), symmetric=True)
    W  = cp.Variable((p,n))

    # region 1 parameters
    L1 = np.array([[2*s]])                      
    M1 = np.array([[1]])

    # region 2 parameters
    L2 = -50000*np.eye(2)
    M2 = np.array([[0, 1],[0, 0]])

    # short form
    S1 = Aa1@Q +Ba@W
    S2 = Aa2@Q +Ba@W

    LMI = [Q >> 0,
            (Aa1@Q +Ba@W +Q@Aa1.T +W.T@Ba.T) << 0,
            (Aa2@Q +Ba@W +Q@Aa2.T +W.T@Ba.T) << 0,
            cp.kron(L1, Q) + cp.kron(M1, S1) + cp.kron(M1.T, S1.T) << 0,
            cp.kron(L1, Q) + cp.kron(M1, S2) + cp.kron(M1.T, S2.T) << 0,
            cp.kron(L2, Q) + cp.kron(M2, S1) + cp.kron(M2.T, S1.T) << -1e-6*np.eye(2*n),
            cp.kron(L2, Q) + cp.kron(M2, S2) + cp.kron(M2.T, S2.T) << -1e-6*np.eye(2*n)]

    pro = cp.Problem(cp.Minimize(0), LMI)
    pro.solve(solver=cp.CVXOPT)
    if pro.status == 'optimal':
        K = (W.value@np.linalg.pinv(Q.value)).flatten()
        print('K =', K)
    else:
        print("Optimization failed:", pro.status)
    return K

def get_tf(Aa1, Aa2, Ba, Ca, Da, Ea, Bq, K, h: int):
    Bk    = Bq*np.vstack((K[1], 0, np.ones((2*h, 1))))
    cl_tf = ss2tf(Aa2+Ba*K, Bk, Ca, Da)
    id_tf = ss2tf(Aa2+Ba*K, Ea, Ca, Da)
    eigv  = np.linalg.eigvals(Aa2+Ba*K)
    print('eig =', eigv)
    return cl_tf, id_tf

def simulate(cl_tf, id_tf, V, R, f, fs=10000, Tf=0.2, st=0.04):
    dt = 1/fs
    t  = np.linspace(0, Tf, int(Tf/dt))

    # response to sine wave input
    r = V*np.sqrt(2)*np.sin(2*math.pi*f*t)
    r[t < st] = 0
    _, y1, _ = signal.lsim((cl_tf.num[0][0], cl_tf.den[0][0]), r, t)

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
    i1[t < st] = 0
    i3[t < st] = 0
    i5[t < st] = 0
    i7[t < st] = 0
    i9[t < st] = 0
    i11[t < st] = 0
    t, y3, x = signal.lsim((id_tf.num[0][0],id_tf.den[0][0]), i3, t)
    t, y5, x = signal.lsim((id_tf.num[0][0],id_tf.den[0][0]), i5, t)
    t, y7, x = signal.lsim((id_tf.num[0][0],id_tf.den[0][0]), i7, t)
    t, y9, x = signal.lsim((id_tf.num[0][0],id_tf.den[0][0]), i9, t)
    t, y11, x = signal.lsim((id_tf.num[0][0],id_tf.den[0][0]), i11, t)

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
    plt.plot(t, r, '--', label='r(t) [V]')
    plt.plot(t, v_o, label='v_o(t) [V]')
    plt.plot(t, i_o, label='i_o(t) [A]')
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.title('Linear load')
    plt.legend()
    plt.grid(True)

    # plot for non-linear load
    plt.figure(2)
    plt.plot(t, r, '--', label='r(t) [V]')
    plt.plot(t, v_o_nl, label='v_o(t) [V]')
    plt.plot(t, i_o_nl, label='i_o(t) [A]')
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.title('Non-linear load')
    plt.legend()
    plt.grid(True)
    plt.show()