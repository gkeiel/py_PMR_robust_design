import PMR_robust_design_functions as ctrl


def main():
    # UPS parameters (LC output filter and load)
    r = 0.2     # filter resistance
    L = 0.001   # filter inductance
    c = 0.00002 # filter capacitance
    R = 6.58    # load resistance
    V = 127     # reference RMS voltage
    f = 50      # reference frequency
    
    # design parameters
    print('Design specifications (ex: h = 3, xi = 0.001, s = 200):')
    h  = int(input('number of resonant modes: '))
    xi = float(input('damping factor in harmonics modes: '))
    s  = int(input('poles minimum real part constraint: '))

    # build plant, controller and augmented models
    A1, A2, B, C, D, E = ctrl.ups_model(r, L, c, R)
    Ar, Br = ctrl.pmr_model(f, h, xi)
    Aa1, Aa2, Bq, Ba, Ca, Da, Ea = ctrl.agm_model(A1, A2, B, C, E, Ar, Br, h)

    # compute augmented state-feedback gain by regional pole placement
    K = ctrl.compute_gain(Aa1, Aa2, Ba, h, s)

    # PMR controller and closed-loop transfer functions
    cl_tf, id_tf = ctrl.get_tf(Aa1, Aa2, Ba, Ca, Da, Ea, Bq, K, h)

    # simulate
    t, r, v_o, i_o, v_o_nl, i_o_nl = ctrl.simulate(cl_tf, id_tf, V, R, f)

    # plot
    ctrl.plot_res(t, r, v_o, i_o, v_o_nl, i_o_nl)


if __name__ == "__main__":
    main()