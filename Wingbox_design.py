import input

wingloading = input.wingloading #import wingloading

b = input.b #import wingspan
S = input.S #import wing area
AR = input.AR #aspect ratio
taper = input.taper #input taper
t_r = input.t_r # input t_r
SMC = b/AR #standar mean chord
widthf = input.widthf

Load_distribution_lift = SMC*wingloading
#load_distribution_wingweight =

print(0.5*(b-widthf), wingloading*S)


