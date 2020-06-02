import input
import matplotlib.pyplot as plt

winglaoding = input.wingloading #import wingloading

b = input.b #import wingspan
S = input.S #import wing area
AR = input.AR #aspect ratio
taper = input.taper #input taper
t_r = input.t_r # input t_r
SMC = b/AR #standar mean chord

lift_d = SMC*winglaoding #not yet adapted to triangle
w_wingd = 1000 # wing weight distribution
w_engine = 10 # engine weight
x_engine = 3 # engine distance from tip

x_loc = []
y_shear = []
y_moment = []

for i in range(1,100):
    x = b*(i/100)
    if x<x_engine:
        Shear_Force = lift_d*x**2/b - w_wingd*x
        Moment = lift_d*x**3/(3*b) - w_wingd*x**2/2
    else:
        Shear_Force = lift_d*x**2/b - w_wingd*x - w_engine
        Moment = lift_d*x**3/(3*b) - w_wingd*x**2/2 - (x-x_engine)*w_engine
    x_loc.append(x)
    y_shear.append(Shear_Force)
    y_moment.append(Moment)


