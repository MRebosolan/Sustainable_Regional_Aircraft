from math import cos, radians
import math
import numpy as np
import matplotlib.pyplot as plt

altitude = 11000
M = 0.8
g = 9.81
range = 2000

H_quarter_sweep = 25

H_MTOW = 36426 * g

H_W_start_cruise = 0.95569551 * H_MTOW
H_W_end_cruise = 0.8587170892 * H_MTOW
H_taper = 0.4
H_S = 81.7  # crj

H_AR = 8
H_t_over_w = 0.4423
H_Swet_Sw = 6.2
H_thrust = H_t_over_w * H_W_start_cruise

H_wingspan = (H_AR * H_S) ** 0.5

# Matteo:
S_quarter_sweep = 25

S_MTOW = 36426 * g

S_W_start_cruise = 0.95569551 * S_MTOW
S_W_end_cruise = 0.8587170892 * S_MTOW
S_taper = 0.4
S_S = 81.7  # crj

S_AR = 12
S_t_over_w = 0.4423
S_Swet_Sw = 5.6
S_thrust = S_t_over_w * S_W_start_cruise

S_wingspan = (S_AR * S_S) ** 0.5

# Tobias:
BWB_quarter_sweep = 45

BWB_MTOW = 35214 * g

BWB_W_start_cruise = BWB_MTOW * 0.95569551
BWB_W_end_cruise = 0.8587170892 * BWB_MTOW
BWB_taper = 0.3
BWB_S = 128.1  # crj

BWB_AR = 4.5  # 5.65
BWB_t_over_w = 0.3892
BWB_Swet_Sw = 3
BWB_thrust = BWB_t_over_w * BWB_W_start_cruise

BWB_wingspan = 26.903

tpaltitude = 7000
tpM = 0.6
g = 9.81
tprange = 2000

tpquarter_sweep = 0

tpMTOW = 36426 * g

tpW_start_cruise = 0.95569551 * tpMTOW
tpW_end_cruise = 0.8587170892 * tpMTOW
tptaper = 0.4
tpS = 81.7  # crj

tpAR = 8
tpt_over_w = 0.4423
tpSwet_Sw = 6.2
tpthrust = tpt_over_w * tpW_start_cruise

tpwingspan = (tpAR * tpS) ** 0.5

crj_quarter_sweep = 26.9

crj_MTOW = 34019 * g

crj_W_start_cruise = 0.95569551 * crj_MTOW
crj_W_end_cruise = 0.8587170892 * crj_MTOW
crj_taper = 0.3
crj_S = 70.6  # crj

crj_AR = 7.6
crj_t_over_w = 1000
crj_Swet_Sw = 6.5
crj_thrust = 61.3 * 1000 * 2

crj_wingspan = (crj_AR * crj_S) ** 0.5


def atmosphere_calculator(h):
    T_grad = -0.0065
    T = 288.15 + T_grad * h
    P = 101325 * (T / 288.15) ** (-9.81 / (T_grad * 287))
    rho = P / (T * 287)
    a = (1.4 * 287 * T) ** 0.5
    return (T, P, rho, a)

def atmosphere_test():
    T, p ,rho ,a = atmosphere_calculator(11000)
    assert T > 216 and T < 217
    assert p > 22614 and p < 22620
    assert rho > 0.36 and rho < 0.37
    assert a > 295 and a < 296
atmosphere_test()


def aerodynamics(quarter_sweep, wingspan, W_start_cruise, W_end_cruise, taper,
                 S, t_over_w, AR, thrust, altitude, M, range, Swet_Sw):


    def CLdes(q, W_start_cruise, W_end_cruise, S, sweep):  # finds design cruise CL for aircraft and airfoil
        CL = 1.1 / q * (0.5 * (W_start_cruise / S + W_end_cruise / S))
        cl = CL / (cos(radians(sweep))) ** 2
        assert CL>0.1 and CL<0.8
        return CL, cl



    def LE_sweep(quarter_sweep, AR, taper):
        return any_sweep(quarter_sweep, AR, taper, 0)

    def any_sweep(quarter_sweep, AR, taper, chord_position):
        quarter_sweep = radians(quarter_sweep)
        any_sweep = np.tan(quarter_sweep) - (4 / AR) * ((0.5 - chord_position) * (1 - taper) / (1 + taper))
        assert any_sweep>0 and any_sweep<90
        return any_sweep

    def oswald_factor(LE_sweep, AR):
        if LE_sweep < 5:
            e = 1.78 * (1 - 0.045 * AR ** 0.68) - 0.64
            assert e>0 and e<1
            return e
        else:
            e = 4.61 * (1 - 0.045 * AR ** 0.68) * (cos(LE_sweep)) ** 0.15 - 3.1
            assert e>0 and e<1
            return e

    def induced_drag(CL, AR=AR, taper=taper, quarter_sweep=quarter_sweep):
        LEsweep = LE_sweep(quarter_sweep, AR, taper)
        e = oswald_factor(LEsweep, AR)
        CD = CL ** 2 / (math.pi * AR * e)
        assert CD>0 and CD<0.2
        return CD

    def CDZERO():
        e = oswald_factor(LE_sweep(quarter_sweep, AR, taper), AR)
        cf = 0.003
        ke = 0.5 * (np.pi * e / cf) ** 0.5
        Emax = ke * (AR / Swet_Sw) ** 0.5
        Cd0 = np.pi * AR * e / (4 * Emax * Emax)
        assert Cd0>0 and Cd0<0.2
        return Cd0

    def drag(V, rho, Cl, S=S, AR=AR):
        Cd0 = CDZERO()
        CD = Cd0 + induced_drag(Cl, AR)
        D = 0.5 * rho * V * V * S * CD
        assert CD>0 and CD<0.2
        return D, CD

    def cruise():
        T, P, rho, a = atmosphere_calculator(altitude)
        V = M * a
        CL = CLdes(0.5 * rho * V * V, W_start_cruise, W_end_cruise, S, quarter_sweep)[0]
        D, CD = drag(V, rho, CL)
        cruise_energy = D * range * 1000  # Joules
        impulse = D * range * 1000 / V
        cruise_power = D * V
        # print("Cruise energy: ", cruise_energy)
        # print("Cruise impulse: ",impulse)

        return cruise_energy, impulse

    def climb(altitude=altitude):
        alt = 0
        t = 0
        impulse = 0
        altlist = []
        thrustlist = []
        v_list = []
        drag_list = []
        climb_list = []
        powerlist = []
        energy = 0
        V_opt = 0

        while alt < altitude and t < 1800:

            start = 80
            speeds = np.arange(start, 220, 1)
            drags = []
            CDS = []
            CLS = []
            for V in speeds:
                T, P, rho, a = atmosphere_calculator(alt)
                CL = 2 * W_start_cruise / (rho * S * V * V)
                D, CD = drag(V, rho, CL)
                drags.append(D)
                CDS.append(CD)
                CLS.append(CL)

            V_previous = V_opt
            D_min = min(drags)
            V_opt = start + drags.index(min(drags))
            CD_opt = CDS[drags.index(min(drags))]
            CL_opt = CLS[drags.index(min(drags))]

            thrust_alt = thrust * (P / 101325) * (288.15 / T) ** 0.5
            thrust_reduced = thrust_alt * .9
            accelerate_power = 0.5 * ((V_opt - V_previous) ** 2) * W_start_cruise / 9.81

            climbrate = ((thrust_reduced - D_min) * V_opt - accelerate_power) / W_start_cruise
            if climbrate < 0:
                climbrate = 0
            climbrate_fpm = climbrate * 60 / 0.3048

            t = t + 1
            alt = alt + climbrate
            impulse = impulse + thrust_reduced

            power = thrust_reduced * V_opt
            powerlist.append(power)
            energy = energy + power
            altlist.append(alt)
            thrustlist.append(thrust_reduced)
            v_list.append(V_opt)
            drag_list.append(D_min)
            climb_list.append(climbrate_fpm)

        time = np.arange(0, len(altlist))

        # plt.plot(time, altlist)
        plt.plot(time, climb_list)
        plt.grid()
        # print("Climb energy: ", energy)
        # print("Climb impulse: ",impulse)

        return energy, impulse



        cr = cruise()
        cli = climb()
        print((cr[1] + cli[1]) / 1000)

        total_energy = cr[0] + cli[0]
        total_impulse = cr[1] + cli[1]

        return total_energy, total_impulse, cr, cli

    print("strutted")
    strutted = aerodynamics(S_quarter_sweep, S_wingspan, S_W_start_cruise, S_W_end_cruise, S_taper,
                            S_S, S_t_over_w, S_AR, S_thrust, altitude, M, range, S_Swet_Sw)

    print("hydrogen")
    hydrogen_jet = aerodynamics(H_quarter_sweep, H_wingspan, H_W_start_cruise, H_W_end_cruise, H_taper,
                                H_S, H_t_over_w, H_AR, H_thrust, altitude, M, range, H_Swet_Sw)
    print("blended")
    blended = aerodynamics(BWB_quarter_sweep, BWB_wingspan, BWB_W_start_cruise, BWB_W_end_cruise, BWB_taper,
                           BWB_S, BWB_t_over_w, BWB_AR, BWB_thrust, altitude, M, range, BWB_Swet_Sw)

    # print(hydrogen[1]/blended[1])
    # print(strutted[1]/blended[1])

    hydrogen_tp = aerodynamics(tpquarter_sweep, tpwingspan, tpW_start_cruise, tpW_end_cruise, tptaper,
                               tpS, tpt_over_w, tpAR, tpthrust, tpaltitude, tpM, tprange, tpSwet_Sw)

    crj = aerodynamics(crj_quarter_sweep, crj_wingspan, crj_W_start_cruise, crj_W_end_cruise, crj_taper,
                       crj_S, crj_t_over_w, crj_AR, crj_thrust, altitude, M, range, crj_Swet_Sw)

    crj_fuel = crj[0] * 290 / 1000 / 1000 / 3600 / range / 75
    print("hydrogen")
    print(hydrogen_jet[0] * 290 / 1000 / 1000 / 3600 / range / 75)

    print("blended")
    print(blended[0] * 290 / 1000 / 1000 / 3600 / range / 75)

    print("strutted")
    print(strutted[0] * 290 / 1000 / 1000 / 3600 / range / 75)

    print("crj")
    print(crj_fuel)

    print("hydrogen")
    print(hydrogen_jet[0] * 290 / 1000 / 1000 / 3600 / range / 75 / crj_fuel)

    print("blended")
    print(blended[0] * 290 / 1000 / 1000 / 3600 / range / 75 / crj_fuel)

    print("strutted")
    print(strutted[0] * 290 / 1000 / 1000 / 3600 / range / 75 / crj_fuel)

