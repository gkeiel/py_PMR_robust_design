import DT_PMR_design_functions as ctrl


def main():
    # UPS parameters (LC output filter and load)
    r = 0.2     # filter resistance
    L = 0.001   # filter inductance
    c = 0.00002 # filter capacitance
    R = 6.58    # nominal load resistance
    V = 127     # reference RMS voltage
    f = 50      # reference frequency
    fs = 20000  # sampling frequency
    T  = 1/fs   # sampling period

    # design parameters
    print('Design specifications (ex: h = 3, xi = 0.001, s = 200):')
    h  = int(input('number of resonant modes: '))
    xi = float(input('damping factor in harmonics modes: '))
    s  = int(input('poles minimum real part constraint: '))

    # build discrete-time (DT) plant, controller and augmented models
    A, B, C, D, E = ctrl.ups_model(r, L, c, R, T)
    Ar, Br = ctrl.pmr_model(f, h, xi, T)
    Aa, Bq, Ba, Ca, Da, Ea = ctrl.agm_model(A, B, C, E, Ar, Br, h)

    # compute DT augmented state-feedback gain
    K = ctrl.compute_gain(Aa, Ba, h, s, T)

    # DT PMR controller and closed-loop transfer functions
    cl_tf, id_tf = ctrl.get_tf(Aa, Ba, Ca, Da, Ea, Bq, K, h, T)

    # simulate
    t, r, v_o, i_o, v_o_nl, i_o_nl = ctrl.simulate(cl_tf, id_tf, V, R, f, T)

    # plot
    ctrl.plot_res(t, r, v_o, i_o, v_o_nl, i_o_nl)


if __name__ == "__main__":
    main()